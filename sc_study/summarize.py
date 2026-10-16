import sys, json
for line in sys.stdin:
    line=line.strip()
    if not line.startswith('{'): 
        print(line); continue
    j=json.loads(line)
    print(f"{j['matrix']:7s} W={j['W']:<6d} off={j['offsets']:7s} dec={j['decoder']:5s} l/d={j['ratio']:.2f} "
          f"succ={j['success']:>3d}/{j['trials']:<3d} silent={j['silent_errors']} undec={j['avg_undecided']:.0f} "
          f"it={j['avg_iters']:.0f} t={j['avg_secs']:.2f}s  H/LB={j['entropy_over_lb']:.3f} CSg/LB={j['cs_global_over_lb']:.3f} CSa/LB={j['cs_adaptive_over_lb']:.3f}")
