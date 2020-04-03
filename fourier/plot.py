import matplotlib.pyplot as plt

with open('fourier_data.txt') as f:
    data = {}
    keys = set()
    correct_ans = {}
    for line in f.readlines():
        line = line.replace('\n','').split(' ')
        a, n, val = int(line[0]), int(line[1]), float(line[2])
        keys.add(a)
        if a not in data:
            data[a] = []
        data[a].append((n, val))
        if a not in correct_ans:
            correct_ans[a] = (n, val)
        else:
            if correct_ans[a][0] < n:
                correct_ans[a] = (n, val)
        #data1.append(int(line[0]))
        #data2.append(float(line[1]))
        #data2.append(abs(float(line[1]) - 0.43434545937892) * pow(10, 10))
    keys = sorted(list(keys), reverse=True)
    def filterData(a):
        data1, data2 = [], []
        for i in data[a]:
            data1.append(i[0])
            #data2.append(i[1])
            data2.append(abs(i[1] - correct_ans[a][1]) * pow(10, 10))
        return data1, data2
    for a in keys:
        data1, data2 = filterData(a)
        plt.plot(data1, data2, label='a = {0}'.format(a / 100.0))
    plt.xscale('log')
    #plt.ylabel('Estimated Probability')
    plt.yscale('log')
    plt.ylabel('Absolute Error $\\times 10^{10}$')
    plt.xlabel('Samples')
    plt.grid()
    plt.legend(loc=4)
    plt.show()