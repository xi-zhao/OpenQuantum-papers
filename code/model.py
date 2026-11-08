"""stories15M (llama2.c 格式) 的加载、分词与 FP32 前向（纯 numpy）。"""
import struct
import numpy as np

MODEL_DIR = "/tmp/claude-0/-home-claude/2b4b55be-e4b8-52fa-a291-d66ebca29c21/scratchpad/model"


class Config:
    def __init__(self, dim, hidden, L, H, Hkv, vocab, seq_len):
        self.dim, self.hidden, self.L, self.H, self.Hkv = dim, hidden, L, H, Hkv
        self.vocab, self.seq_len = abs(vocab), seq_len
        self.shared = vocab > 0
        self.hs = dim // H


def load(path=f"{MODEL_DIR}/stories15M.bin"):
    raw = open(path, "rb").read()
    cfg = Config(*struct.unpack("7i", raw[:28]))
    f = np.frombuffer(raw, dtype=np.float32, offset=28)
    d, hd, L, V = cfg.dim, cfg.hidden, cfg.L, cfg.vocab
    kv = cfg.Hkv * cfg.hs
    o = 0

    def take(*shape):
        nonlocal o
        n = int(np.prod(shape))
        a = f[o:o + n].reshape(shape).astype(np.float64)
        o += n
        return a

    W = {}
    W["tok"] = take(V, d)
    W["rms_att"] = take(L, d)
    W["wq"] = take(L, d, d)
    W["wk"] = take(L, kv, d)
    W["wv"] = take(L, kv, d)
    W["wo"] = take(L, d, d)
    W["rms_ffn"] = take(L, d)
    W["w1"] = take(L, hd, d)
    W["w2"] = take(L, d, hd)
    W["w3"] = take(L, hd, d)
    W["rms_final"] = take(d)
    take(cfg.seq_len, cfg.hs // 2)  # freq_cis_real（旧格式，未使用）
    take(cfg.seq_len, cfg.hs // 2)  # freq_cis_imag
    W["cls"] = W["tok"] if cfg.shared else take(V, d)
    return cfg, W


def rope_tables(cfg, T):
    """与 llama2.c run.c 一致：按维度对 (i, i+1) 旋转，频率由 i % head_size 决定。"""
    i = np.arange(0, cfg.dim, 2)
    head_dim = i % cfg.hs
    freq = 1.0 / (10000.0 ** (head_dim / cfg.hs))
    ang = np.arange(T)[:, None] * freq[None, :]
    return np.cos(ang), np.sin(ang)          # (T, dim/2)


def apply_rope(x, cos, sin):
    x0, x1 = x[:, 0::2], x[:, 1::2]
    out = np.empty_like(x)
    out[:, 0::2] = x0 * cos - x1 * sin
    out[:, 1::2] = x0 * sin + x1 * cos
    return out


def rmsnorm(x, w, eps=1e-5):
    return x / np.sqrt((x * x).mean(-1, keepdims=True) + eps) * w


def forward_fp(cfg, W, tokens, record=None):
    """整段序列的因果前向，返回 (T, vocab) logits。record: dict，可收集中间激活用于校准。"""
    T = len(tokens)
    x = W["tok"][tokens].copy()
    cos, sin = rope_tables(cfg, T)
    mask = np.triu(np.full((T, T), -np.inf), 1)
    for l in range(cfg.L):
        xn = rmsnorm(x, W["rms_att"][l])
        q, k, v = xn @ W["wq"][l].T, xn @ W["wk"][l].T, xn @ W["wv"][l].T
        q, k = apply_rope(q, cos, sin), apply_rope(k, cos, sin)
        att_out = np.zeros_like(q)
        for h in range(cfg.H):
            sl = slice(h * cfg.hs, (h + 1) * cfg.hs)
            s = q[:, sl] @ k[:, sl].T / np.sqrt(cfg.hs) + mask
            s = np.exp(s - s.max(-1, keepdims=True))
            p = s / s.sum(-1, keepdims=True)
            att_out[:, sl] = p @ v[:, sl]
        x = x + att_out @ W["wo"][l].T
        xn2 = rmsnorm(x, W["rms_ffn"][l])
        h1, h3 = xn2 @ W["w1"][l].T, xn2 @ W["w3"][l].T
        g = h1 / (1 + np.exp(-h1)) * h3
        x = x + g @ W["w2"][l].T
        if record is not None:
            for name, val in (("xn", xn / W["rms_att"][l]), ("q", q), ("k", k), ("v", v), ("att", att_out),
                              ("xn2", xn2 / W["rms_ffn"][l]), ("h1", h1), ("h3", h3), ("g", g), ("res", x)):
                record.setdefault((name, l), []).append(val)
    xf = rmsnorm(x, W["rms_final"])
    if record is not None:
        record.setdefault(("xf", 0), []).append(xf / W["rms_final"])
    return xf @ W["cls"].T


# ---------------- 分词器（llama2.c tokenizer.bin） ----------------
class Tokenizer:
    def __init__(self, path=f"{MODEL_DIR}/tokenizer.bin", vocab=32000):
        b = open(path, "rb").read()
        self.max_len = struct.unpack("i", b[:4])[0]
        o, self.vocab, self.scores = 4, [], []
        for _ in range(vocab):
            sc, n = struct.unpack("fi", b[o:o + 8]); o += 8
            self.vocab.append(b[o:o + n]); o += n
            self.scores.append(sc)
        self.lookup = {t: i for i, t in enumerate(self.vocab)}

    def encode(self, text, bos=True):
        ids = []
        if text:
            ids.append(self.lookup[b" "])
        for ch in text:
            bs = ch.encode("utf-8")
            if bs in self.lookup:
                ids.append(self.lookup[bs])
            else:
                ids.extend(c + 3 for c in bs)
        while True:
            best, bi, bid = -1e10, -1, -1
            for i in range(len(ids) - 1):
                m = self.vocab[ids[i]] + self.vocab[ids[i + 1]]
                j = self.lookup.get(m)
                if j is not None and self.scores[j] > best:
                    best, bi, bid = self.scores[j], i, j
            if bi < 0:
                break
            ids[bi:bi + 2] = [bid]
        return ([1] if bos else []) + ids

    def decode(self, ids):
        out = b""
        prev = None
        for t in ids:
            if t in (1, 2):
                prev = t
                continue
            p = self.vocab[t]
            if prev == 1 and p.startswith(b" "):
                p = p[1:]
            if p.startswith(b"<0x") and p.endswith(b">") and len(p) == 6:
                p = bytes([int(p[3:5], 16)])
            out += p
            prev = t
        return out.decode("utf-8", errors="replace")
