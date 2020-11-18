def suffixArray(s):
    satups = sorted([(s[i:],i) for i in range(0, len(s))])
    print(satups)
    return map(lambda x:  x[1],satups)

def bwtViaSa(t):
    bw =[]
    for si in suffixArray(t):
        if si ==0: bw.append("$")
        else:bw.append(t[si-1])
    return "".join(bw)

print(bwtViaSa('SAAASSSSASG'))