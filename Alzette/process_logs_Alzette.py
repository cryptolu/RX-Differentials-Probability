from process_logs import LogFile

TABLE_PREFIX = r"""
\begin{table}[htbp]
\begin{center}
\begin{tabular}{|c|c|c|c||c|c|c|}
\hline
& \multicolumn{3}{c||}{NEQ optimization} & \multicolumn{3}{c|}{Probability optimization} \\ \hline
$k$ & nb NEQ & best LB & wt & nb NEQ &  best LB & wt \\ \hline
%%"""

TABLE_SUFFIX = r"""%%
\hline
\end{tabular}
\caption{
\TabAlzPref{%(hours)dh} $c_%(const_id)d =$ \texttt{0x%(const)08x}.
\TabAlzSuf{} %(has_no_match)s}
\label{tab:alzettec%(const_id)d}
\end{center}
\end{table}
"""


CONSTS = [
    0xb7e15162,
    0xbf715880,
    0x38b4da56,
    0x324e7738,
    0xbb1185eb,
    0x4f7c7b57,
    0xcfbfa1c8,
    0xc2b3293d,
]

best_k1 = [(999, 0)] * 8
best_k2plus = [(999, 0)] * 8

for CONST_ID in range(8):
    TL_NEQ = 3600
    if CONST_ID == 0:
        TL_FULL = 14400
    else:
        TL_FULL = 7200

    HOURS = (TL_FULL + 3600) // 3600

    CONST = CONSTS[CONST_ID]

    best_prob_neq = 999
    best_prob_full = 999
    best_neq_neq = 999
    best_neq_full = 999
    datas = []
    for k in range(1, 17):
        log_fname = f"logs/alzette_const{CONST:08x}_rotation{k:d}_neqopt_time{TL_NEQ}.txt"
        trail_fname = f"trails/alzette_const{CONST:08x}_rotation{k:d}_neqopt_time{TL_NEQ}.txt"
        FN = LogFile(log_fname)
        with open(trail_fname) as f:
            neq, prob_neq = f.read().splitlines()[2].strip().split()
            neq = int(neq)
            prob_neq = -float(prob_neq)


        log_fname = f"logs/alzette_const{CONST:08x}_rotation{k:d}_fullopt_time{TL_FULL}.txt"
        trail_fname = f"trails/alzette_const{CONST:08x}_rotation{k:d}_fullopt_time{TL_FULL}.txt"
        FF = LogFile(log_fname)
        with open(trail_fname) as f:
            neq_full, prob_full = f.read().splitlines()[2].strip().split()
            neq_full = int(neq_full)
            prob_full = -float(prob_full)

        lb_neq, num_neq = FN.lowerX[-1], FN.upperY[-1]
        lb_full, num_full = FF.lowerX[-1], FF.upperY[-1]
        assert num_neq == neq
        assert abs(num_full - prob_full) < 1e-3, (num_full, prob_full)

        has_no_match_neq = FN.has_no_match
        has_no_match_full = FF.has_no_match
        assert has_no_match_full is not None
        assert has_no_match_neq is None

        datas.append((k, neq, prob_neq, neq_full, prob_full, lb_neq, num_neq, lb_full, num_full, has_no_match_neq, has_no_match_full))

        best_prob_neq = min(best_prob_neq, prob_neq)
        best_prob_full = min(best_prob_full, num_full)

        best_neq_neq = min(best_neq_neq, neq)
        best_neq_full = min(best_neq_full, neq_full)

        assert best_prob_full - 1e-3 < best_prob_neq
        time_used = FN.timesY[-1] + FF.timesY[-1]
        if k == 1:
            best_k1[CONST_ID] = min(best_k1[CONST_ID], (num_full, time_used))
        elif k > 1:
            best_k2plus[CONST_ID] = min(best_k2plus[CONST_ID], (num_full, time_used, k))

    has_no_match = False
    print(TABLE_PREFIX % dict(const=CONST, const_id=CONST_ID))
    for data in datas:
        k, neq, prob_neq, neq_full, prob_full, lb_neq, num_neq, lb_full, num_full, has_no_match_neq, has_no_match_full = data
        # has_no_match |= has_no_match_neq
        has_no_match |= has_no_match_full
        if abs(lb_neq - num_neq) < 0.1:
            print("OPTIMAL NEQ!!!")
            raise
        if abs(lb_full - num_full) < 0.01:
            print("OPTIMAL FULL!!!")
            raise

        #print(best_prob_neq, prob_neq, abs(prob_neq - best_prob_neq))
        pref1 = r"\underline{" if abs(prob_neq - best_prob_neq) < 0.01 else "{"
        pref2 = r"\underline{" if abs(prob_full - best_prob_full) < 0.01 else "{"
        pref3 = r"\underline{" if abs(neq - best_neq_neq) < 0.01 else "{"
        pref4 = r"\underline{" if abs(neq_full - best_neq_full) < 0.01 else "{"
        suff1 = "" # r"\textdagger" if FN.has_no_match else ""
        suff2 = r"\textdagger" if has_no_match_full else ""
        print(rf"{k} & {pref3}{num_neq:.0f}}} & {lb_neq:.0f} & {pref1}{prob_neq:.2f}}}{suff1} & {pref4}{neq_full}}} & {lb_full:.2f} & {pref2}{num_full:.2f}}}{suff2} \\")
    print("%")

    has_no_match = r"\HasNoMatch{}" if has_no_match else ""

    print(TABLE_SUFFIX % dict(const=CONST, const_id=CONST_ID, hours=HOURS, has_no_match=has_no_match))

print("------")
for ci in range(8):
    row = best_k1[ci]
    #print(row)
    print("c%d" % ci)
    print(f"& {row[0]:.2f} & {row[1]/60.0:.2f} min", end=" ")
    row = best_k2plus[ci]
    #print()
    #print(row)
    print(f"& {row[0]:.2f} & {row[1]/60.0:.2f} min")
    print()

