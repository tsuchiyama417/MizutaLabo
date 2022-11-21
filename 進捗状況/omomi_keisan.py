import linecache
path_aaindex1 = 'C:/Users/k/Dropbox/kenkyu/M1/database/d/aaindex1.txt'
def get_value(s):
    s = 'H ' + s
    flag = 0
    row = 0
    with open(path_aaindex1, 'r') as f:
        for line in f:
            if line == s:
                flag = 1
            if ("I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V" in line) and flag == 1:
                return row + 2
            row += 1

xyz = [[0,0,0] for _ in range(20)]

cnt = 0
with open('coordinate.txt',  'r') as f:
    for line in f:
        r = get_value(line)
        l = []
        l.append(linecache.getline(path_aaindex1, r).split())
        l.append(linecache.getline(path_aaindex1, r + 1).split())

        xyz[0][cnt] = float(l[0][0]) # A
        xyz[1][cnt] = float(l[0][4]) # C
        xyz[2][cnt] = float(l[0][3]) # D
        xyz[3][cnt] = float(l[0][6]) # E
        xyz[4][cnt] = float(l[1][3]) # F
        xyz[5][cnt] = float(l[0][7]) # G
        xyz[6][cnt] = float(l[0][8]) # H
        xyz[7][cnt] = float(l[0][9]) # I
        xyz[8][cnt] = float(l[1][1]) # K
        xyz[9][cnt] = float(l[1][0]) # L
        xyz[10][cnt] = float(l[1][2]) # M
        xyz[11][cnt] = float(l[0][2]) # N
        xyz[12][cnt] = float(l[1][4]) # P
        xyz[13][cnt] = float(l[0][5]) # Q
        xyz[14][cnt] = float(l[0][1]) # R
        xyz[15][cnt] = float(l[1][5]) # S
        xyz[16][cnt] = float(l[1][6]) # T
        xyz[17][cnt] = float(l[1][9]) # V
        xyz[18][cnt] = float(l[1][7]) # W
        xyz[19][cnt] = float(l[1][8]) # Y
        
        cnt += 1

#アルファベット順
omomilist = [1.083, 1.859, 1.263, 1.172, 1.413, 1.150, 1.643, 1.228, 1.236, 1.015,
            1.617, 1.391, 1.324, 1.405, 1.257, 1.178, 1.271, 1.163, 1.959, 1.535]

for i in range(20):
    for j in range(3):
            xyz[i][j] = xyz[i][j] * omomilist[i]

dic = {"A":(xyz[0]), "C":(xyz[1]), "D":(xyz[2]), "E":(xyz[3]), "F":(xyz[4]),
        "G":(xyz[5]), "H":(xyz[6]), "I":(xyz[7]), "K":(xyz[8]), "L":(xyz[9]),
        "M":(xyz[10]), "N":(xyz[11]), "P":(xyz[12]), "Q":(xyz[13]), "R":(xyz[14]),
        "S":(xyz[15]), "T":(xyz[16]), "V":(xyz[17]), "W":(xyz[18]), "Y":(xyz[19])}