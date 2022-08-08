import numpy as np
import matplotlib.pyplot as plt

animal = ("Human","Gorilla","Pygmi_chimpanzee","Common_chimpanzee","Fin_whale","Blue_whale","Rat","Mouse","Opossum")

def pl(Arr):
    
    # グラフの枠を作成
    fig = plt.figure()
    #色の取得
    colors = ["red","coral","gold","olive","cyan","blue","purple","magenta","green"]

    # X_listには各生物種のx座標が格納されている。Y,Zも同様。
    xyzArr = Arr
    X_list = []
    Y_list = []
    Z_list = []
    for i in range(0,27,3):
        X_list.append(xyzArr[i])
        Y_list.append(xyzArr[i+1])
        Z_list.append(xyzArr[i+2])

    # S = 分散共分散行列
    S_list = []
    for i in range(9):
        Sxx = np.var(X_list[i])
        Syy = np.var(Y_list[i])
        Szz = np.var(Z_list[i])
        Sxy = np.cov(X_list[i], Y_list[i])
        Sxz = np.cov(X_list[i], Z_list[i])
        Syz = np.cov(Y_list[i], Z_list[i])
        S = np.array([[Sxx,      Sxy[0][1], Sxz[0][1]],
                     [Sxy[0][1], Syy,       Syz[0][1]],
                     [Sxz[0][1], Syz[0][1], Szz]])
        S_list.append(S)

  
    # 固有値・固有ベクトルを計算
    vector = []
    eigenvalue = []
    count = 0
    for i in S_list:
        ans = np.linalg.eig(i)
        print(animal[count])
        print("分散共分散行列")
        print(i)
        print("固有値")
        print(ans[0]) # 固有値
        print("固有ベクトル")
        # 固有ベクトルのz座標が負の値のとき固有ベクトルの向きを反転する
        if ans[1][2][0] < 0:
            ans[1][0][0] = ans[1][0][0] * -1
            ans[1][1][0] = ans[1][1][0] * -1
            ans[1][2][0] = ans[1][2][0] * -1
        elif ans[1][2][1] < 0:
            ans[1][0][1] = ans[1][0][1] * -1
            ans[1][1][1] = ans[1][1][1] * -1
            ans[1][2][1] = ans[1][2][1] * -1
        elif ans[1][2][2] < 0:
            ans[1][0][2] = ans[1][0][2] * -1
            ans[1][1][2] = ans[1][1][2] * -1
            ans[1][2][2] = ans[1][2][2] * -1
        
        print(ans[1]) # 固有ベクトル
     
        vector.append(ans[1])
        eigenvalue.append(ans[0])
        count += 1

    # 3次元の点を2次元の点へ変換
    x = []
    y = []
    
    ### 生物種に応じて変更 ###
    i = 0
    ####                  ####

    # 第一主成分、第二主成分、第三主成分の大きさを比較する必要あり
    # 各自コードを改変してください
    for j in range(len(X_list[i])):
        # 第一主成分
        x.append(vector[i][0][0] * X_list[i][j] + vector[i][1][0] * Y_list[i][j] + vector[i][2][0] * Z_list[i][j]) 
        # 第二主成分
        if eigenvalue[i][1] > eigenvalue[i][2]:
            y.append(vector[i][0][1] * X_list[i][j] + vector[i][1][1] * Y_list[i][j] + vector[i][2][1] * Z_list[i][j])
        else:
            y.append(vector[i][0][2] * X_list[i][j] + vector[i][1][2] * Y_list[i][j] + vector[i][2][2] * Z_list[i][j])

    for k in range(len(X_list[i])):
        plt.scatter(x[k], y[k], color = colors[i])

    # タイトルの設定
    fig.suptitle(animal[i], size=18, weight=2)
    # 最後に.show()を書いてグラフ表示
    plt.show()



