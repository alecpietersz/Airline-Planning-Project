import itertools

hub = 'a'

s = [hub, hub, 'b', 'c', 'd']
res = set()
for k in range(1, len(s)+1):
    for comb in itertools.permutations(s, k):
            if comb[0] == hub and comb[-1] == hub and len(comb)>2:
                print(comb)
                res.add(''.join(comb))

print(res)