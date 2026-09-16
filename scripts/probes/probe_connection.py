"""Read-only ProxyJump connectivity check: laptop -> jump host -> target compute host.

Opens a direct-tcpip tunnel through the jump host and authenticates to the
target over it, then runs a few harmless read-only commands. Touches nothing.
All connection details come from config.txt (gitignored); nothing hardcoded.
"""
import os
import sys
import paramiko

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def _load_config():
    cfg = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "config.txt")
    vals = {}
    if os.path.exists(cfg):
        with open(cfg, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if ":" in line:
                    k, v = line.split(":", 1)
                    vals[k.strip().lower()] = v.strip()
    return vals

_cfg = _load_config()
if not _cfg:
    raise SystemExit(
        "config.txt not found. Copy config.example.txt to config.txt (project root) "
        "and fill in host/user/pwd/atq/tgt_host (gitignored, never commit it)."
    )

JUMP_HOST = _cfg["host"]
JUMP_USER = _cfg["user"]
JUMP_PASS = _cfg["pwd"]

TGT_HOST = _cfg["tgt_host"]
TGT_PORT = 22
TGT_USER = _cfg["user"]
TGT_PASS = _cfg["atq"]

RO_CMD = "hostname; whoami; echo '--- uname ---'; uname -a; echo '--- uptime ---'; uptime; echo '--- pwd ---'; pwd"


def main():
    jump = paramiko.SSHClient()
    jump.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        jump.connect(JUMP_HOST, username=JUMP_USER, password=JUMP_PASS,
                     timeout=25, allow_agent=False, look_for_keys=False)
        print("[1/2] jump OK  -> %s (%s)" % (JUMP_HOST, JUMP_USER))
    except Exception as ex:
        print("[1/2] jump FAILED -> %s: %s" % (type(ex).__name__, ex))
        return

    try:
        transport = jump.get_transport()
        chan = transport.open_channel(
            "direct-tcpip", (TGT_HOST, TGT_PORT), ("127.0.0.1", 0), timeout=25)
    except Exception as ex:
        print("[tunnel] FAILED to open channel to %s:%d -> %s: %s"
              % (TGT_HOST, TGT_PORT, type(ex).__name__, ex))
        jump.close()
        return

    tgt = paramiko.SSHClient()
    tgt.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        tgt.connect(TGT_HOST, port=TGT_PORT, username=TGT_USER, password=TGT_PASS,
                    sock=chan, timeout=25, allow_agent=False, look_for_keys=False)
        print("[2/2] quantumatk hop OK -> %s (%s)\n" % (TGT_HOST, TGT_USER))
    except paramiko.AuthenticationException as ex:
        print("[2/2] quantumatk AUTH FAILED -> %s (wrong password for this hop?)" % ex)
        jump.close()
        return
    except Exception as ex:
        print("[2/2] quantumatk connect FAILED -> %s: %s" % (type(ex).__name__, ex))
        jump.close()
        return

    try:
        _, out, err = tgt.exec_command(RO_CMD, timeout=30)
        print(out.read().decode(errors="replace"))
        e = err.read().decode(errors="replace")
        if e.strip():
            print("STDERR:", e)
    finally:
        tgt.close()
        jump.close()
    print("=== connection closed cleanly, nothing modified ===")


if __name__ == "__main__":
    main()
