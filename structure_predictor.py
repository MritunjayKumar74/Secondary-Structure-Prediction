P_alpha = {
    'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20,
    'Q': 1.17, 'W': 1.14, 'V': 1.14, 'F': 1.12, 'K': 1.07,
    'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79,
    'C': 0.77, 'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53
}

P_beta = {
    'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29,
    'F': 1.28, 'Q': 1.23, 'L': 1.22, 'T': 1.20, 'W': 1.19,
    'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80, 'K': 0.74,
    'S': 0.72, 'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26
}

possible_sequences = []

sequence = ("MAQWNQLQQLDTRYLEQLHQLYSDSFPMELRQFLAPWIESQDWAYAASKESHATLVFHNLLGE"
            "IDQQYSRFLQESNVLYQHNLRRIKQFLQSRYLEKPMEIARIVARCLWEESRLLQTAATAAQQG"
            "GQANHPTAAVVTEKQQMLEQHLQDVRKRVQDLEQKMKVVENLQDDFDFNYKTLKSQGDMQDLN"
            "GNNQSVTRQKMQQLEQMLTALDQMRRSIVSELAGLLSAMEYVQKTLTDEELADWKRRQQIACI"
            "GGPPNICLDRLENWITSLAESQLQTRQQIKKLEELQQKVSYKGDPIVQHRPMLEERIVELFRN"
            "LMKSAFVVERQPCMPMHPDRPLVIKTGVQFTTKVRLLVKFPELNYQLKIKVCIDKDSGDVAAL"
            "RGSRKFNILGTNTKVMNMEESNNGSLSAEFKHLTLREQRCGNGGRANCDASLIVTEELHLITF"
            "ETEVYHQGLKIDLETHSLPVVVISNICQMPNAWASILWYNMLTNNPKNVNFFTKPPIGTWDQV"
            "AEVLSWQFSSTTKRGLSIEQLTTLAEKLLGPGVNYSGCQITWAKFCKENMAGKGFSFWVWLDN"
            "IIDLVKKYILALWNEGYIMGFISKERERAILSTKPPGTFLLRFSESSKEGGVTFTWVEKDISG"
            "KTQIQSVEPYTKQQLNNMSFAEIIMGYKIMDATNILVSPLVYLYPDIPKEEAFGKYCRPESQE"
            "HPEADPGSAAPYLKTKFICVTPTTCSNTIDLPMSPRTLDSLMQFGNNGEGAEPSAGGQFESLT"
            "FDMELTSECATSPM")

def find_helix_nucleation_sites(sequence):
    nucleation_sites = []

    for i in range(len(sequence) - 5):
        part = sequence[i:i+6]

        count = 0

        for x in part:
            if P_alpha[x] > 1:
                count+=1

        if count >= 4:
            nucleation_sites.append((i, i+6))

    return nucleation_sites

def extend_helix(sequence, start, end):
    while end < len(sequence):

        part = sequence[end-3: end+1]

        if len(part) < 4:
            break

        sum = 0

        for x in part:
            sum += P_alpha[x]

        if sum >= 4.0:
            end += 1
        else:
            break

    while start > 0:

        part = sequence[start-1: start+3]

        if len(part) < 4:
            break

        sum = 0

        for x in part:
            sum += P_alpha[x]

        if sum >= 4.0:
            start -= 1
        else:
            break

    return start, end

def merge_overlaps(helix_set):

    n = len(helix_set)

    result = []

    for i in range(n):
        x = helix_set[i][0]
        y = helix_set[i][1]

        if result and result[-1][1] >= y:
            continue

        for j in range(i+1, n):
            if helix_set[j][0] < y:
                y = max(y, helix_set[j][1])


        result.append((x, y))

    return result

def print_helix(helix_set, sequence):
    print()
    print("********************Part (a)********************")
    print()
    print("The sequence regions that are helical in nature: ")
    print()

    for x, y in helix_set:
        print(sequence[x: y] + " (Position " + str(x) + "-" + str(y) + ")")

    print()
    print("********************Part (a) end********************")
    print()

def find_sheet_nucleation_sites(sequence):
    nucleation_sites = []

    for i in range(len(sequence) - 4):
        part = sequence[i:i+5]

        count = 0

        for x in part:
            if P_beta[x] > 1:
                count+=1

        if count >= 3:
            nucleation_sites.append((i, i+5))

    return nucleation_sites

def extend_sheet(sequence, start, end):
    while end < len(sequence):

        part = sequence[end-3: end+1]

        if len(part) < 4:
            break

        sum = 0

        for x in part:
            sum += P_beta[x]

        if sum >= 4.0:
            end += 1
        else:
            break

    while start > 0:

        part = sequence[start-1: start+3]

        if len(part) < 4:
            break

        sum = 0

        for x in part:
            sum += P_beta[x]

        if sum >= 4.0:
            start -= 1
        else:
            break

    return start, end

def print_sheet(sheet_set, sequence):
    print()
    print("********************Part (b)********************")
    print()
    print("The sequence regions that are sheet in nature: ")
    print()

    for x, y in sheet_set:
        print(sequence[x: y] + " (Position " + str(x) + "-" + str(y) + ")")

    print()
    print("********************Part (b) end********************")
    print()

def find_conflicts(helix_set, sheet_set):
    conflicts = []
    for h_start, h_end in helix_set:
        for s_start, s_end in sheet_set:
            overlap_start = max(h_start, s_start)
            overlap_end = min(h_end, s_end)
            if overlap_start < overlap_end:
                conflicts.append((overlap_start, overlap_end))
    return conflicts

def print_conflict(conflict_set, sequence):
    print()
    print("********************Part (c)********************")
    print()
    print("The sequence regions that are conflicting: ")
    print()

    for x, y in conflict_set:
        print(sequence[x: y] + " (Position " + str(x) + "-" + str(y) + ")")

    print()
    print("********************Part (c) end********************")
    print()

def solve_conflict(helix_set, sheet_set, conflict_set, sequence):

    final_structure = ["-"]*len(sequence)

    final_value = {}

    for start, end in conflict_set:
        alpha_sum = 0
        beta_sum = 0

        for x in sequence[start: end]:
            alpha_sum += P_alpha[x]
            beta_sum += P_beta[x]

        if alpha_sum > beta_sum:
            final_value[(start, end)] = "H"
        else:
            final_value[(start, end)] = "E"

    for start, end in helix_set:
        for i in range(start, end):
            final_structure[i] = "H"

    for start, end in sheet_set:
        for i in range(start, end):
            final_structure[i] = "E"

    for (start, end), value in final_value.items():
        for i in range(start, end):
            final_structure[i] = value

    return final_structure, final_value

helix_set = set()

for x, y in find_helix_nucleation_sites(sequence):
    helix_set.add(extend_helix(sequence, x, y))

helix_set = list(sorted(helix_set))
helix_set = merge_overlaps(helix_set)

print_helix(helix_set, sequence)

sheet_set = set()

for x, y in find_sheet_nucleation_sites(sequence):
    sheet_set.add(extend_sheet(sequence, x, y))

sheet_set = list(sorted(sheet_set))
sheet_set = merge_overlaps(sheet_set)

print_sheet(sheet_set, sequence)

conflict_set = find_conflicts(helix_set, sheet_set)

print_conflict(conflict_set, sequence)

final_structure, final_value = solve_conflict(helix_set, sheet_set, conflict_set, sequence)

structure_string = ''.join(final_structure)

print("Final complete sequence: ")
print()

for i in range(0, len(sequence), 60):
    print(sequence[i:i+60])
    print(structure_string[i:i+60])
    print()


# I have used these rules -

# For identifying nucleation sites:
# Helix: 4 out of 6 consecutive residues with P(H) > 1
# Sheet: 3 out of 5 consecutive residues with P(S) > 1

# For extending:
# Extend in both directions using a 4-residue window
# Continue extending if the sum of P(H) or P(S) in the window ≥ 4.0
# Stop when the sum < 4.0
