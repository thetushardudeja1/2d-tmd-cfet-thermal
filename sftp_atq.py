"""SFTP push/pull for a compute host reached through an SSH jump host (ProxyJump).

Same ProxyJump tunnel as ssh_run_atq.py. Read/write of files only; it never
deletes or renames anything on the remote host. An upload whose basename
collides with an existing remote file is REFUSED unless --force is given.

    python sftp_atq.py put <remote_dir> <local_path> [<local_path> ...] [--force]
    python sftp_atq.py get <local_dir>  <remote_path> [<remote_path> ...]

Remote paths may be absolute or ~-relative (expanded from config.txt's user).
All connection details are read from ./config.txt (gitignored); nothing here
is hardcoded.
"""
import os
import posixpath
import sys
import paramiko

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# --- filled in by load_config() from config.txt (gitignored, not in this repo) ---
JUMP_HOST = None
JUMP_USER = None
JUMP_PASS = None
TGT_HOST = None
TGT_PORT = 22
TGT_USER = None
TGT_PASS = None
HOME = None


def load_config():
    global JUMP_HOST, JUMP_USER, JUMP_PASS, TGT_HOST, TGT_USER, TGT_PASS, HOME
    cfg = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.txt")
    if not os.path.exists(cfg):
        raise SystemExit(
            "config.txt not found. Copy config.example.txt to config.txt and "
            "fill in host/user/pwd/atq/tgt_host/home (this file is gitignored, never commit it)."
        )
    vals = {}
    with open(cfg, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if ":" in line:
                k, v = line.split(":", 1)
                vals[k.strip().lower()] = v.strip()
    JUMP_HOST = vals.get("host", JUMP_HOST)
    JUMP_USER = vals.get("user", JUMP_USER)
    JUMP_PASS = vals.get("pwd", JUMP_PASS)
    TGT_HOST = vals["tgt_host"]
    TGT_USER = vals.get("user", JUMP_USER)
    TGT_PASS = vals.get("atq", TGT_PASS)
    HOME = vals.get("home", "/home/%s" % JUMP_USER)


def open_sftp():
    jump = paramiko.SSHClient()
    jump.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    jump.connect(JUMP_HOST, username=JUMP_USER, password=JUMP_PASS,
                 timeout=25, allow_agent=False, look_for_keys=False)
    chan = jump.get_transport().open_channel(
        "direct-tcpip", (TGT_HOST, TGT_PORT), ("127.0.0.1", 0), timeout=25)
    tgt = paramiko.SSHClient()
    tgt.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    tgt.connect(TGT_HOST, port=TGT_PORT, username=TGT_USER, password=TGT_PASS,
                sock=chan, timeout=25, allow_agent=False, look_for_keys=False)
    return jump, tgt, tgt.open_sftp()


def expand(p):
    return p.replace("~", HOME, 1) if p.startswith("~") else p


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    load_config()
    mode = sys.argv[1]
    force = "--force" in sys.argv
    args = [a for a in sys.argv[2:] if a != "--force"]

    jump, tgt, sftp = open_sftp()
    try:
        if mode == "put":
            remote_dir, locals_ = expand(args[0]), args[1:]
            try:
                sftp.stat(remote_dir)
            except IOError:
                sftp.mkdir(remote_dir)
                print("mkdir %s" % remote_dir)
            for lp in locals_:
                dst = posixpath.join(remote_dir, os.path.basename(lp))
                exists = True
                try:
                    sftp.stat(dst)
                except IOError:
                    exists = False
                if exists and not force:
                    print("REFUSED (exists, use --force): %s" % dst)
                    continue
                sftp.put(lp, dst)
                print("%10d  %s" % (os.path.getsize(lp), dst))
        elif mode == "get":
            local_dir, remotes = args[0], args[1:]
            os.makedirs(local_dir, exist_ok=True)
            for r in remotes:
                r = expand(r)
                dst = os.path.join(local_dir, posixpath.basename(r))
                sftp.get(r, dst)
                print("%10d  %s" % (os.path.getsize(dst), dst))
        else:
            print("first arg must be 'put' or 'get'")
    finally:
        sftp.close()
        tgt.close()
        jump.close()


if __name__ == "__main__":
    main()
