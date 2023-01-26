def trypsin(sequence, missed_clevage, min_len, max_len):
	if 'K' in sequence or 'R' in sequence:
		get_dup_k = [i for i in range(len(sequence)) if sequence.startswith('K', i)]
		get_dup_r = [j for j in range(len(sequence)) if sequence.startswith('R', j)]
		merge_list = sorted(get_dup_k + get_dup_r)
		merge_list_fltrd = [i for i in merge_list if i+1 < len(sequence) and sequence[i + 1] !='P'] #look for KP or RP position
		merge_list_fltrd.append(len(sequence))
		initialize = 0
		for iter_lst in range(len(merge_list_fltrd) - int(missed_clevage)):
			peptide = (sequence[initialize: int(merge_list_fltrd[iter_lst + missed_clevage]) + 1])
			if len(peptide) >= int(min_len) and len(peptide) <= int(max_len):
				yield peptide
			initialize = merge_list_fltrd[iter_lst] + 1

'''
a = trypsin("MGKNKLLHPSLVLLLLVLLPTDASVSGKPQYMVLVPSLLHTETTEKGCVLLSYLNETVTVSASLESVRGNRSLFTDLEAENDVLHCVAFAVPKSSSNEEVMFLTVQVKGPTQEFKKRTTVMVKNEDSLVFVQTDKSIYKPGQTVKFRVVSMDENFHPLNELIPLVYIQDPKGNRIAQWQSFQLEG", 0, 7, 25)
for i in a:
    print (i)
'''
