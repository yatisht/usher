import matplotlib.pyplot as plt
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot recombinants')
    parser.add_argument("-v", type=str,
                        help="VCF containing donor, acceptor and recombinant"
                        " node")
    parser.add_argument("-l", type=int,
                        help="genome length")
    parser.add_argument("-s1", type=int,
                        help="start low")
    parser.add_argument("-s2", type=int,
                        help="start high")
    parser.add_argument("-e1", type=int,
                        help="end low")
    parser.add_argument("-e2", type=int,
                        help="end high")


    args = vars(parser.parse_args())
    vcf_filename = args.get('v', '')
    length = args.get('l', '')
    s1 = args.get('s1', '')
    e1 = args.get('e1', '')
    s2 = args.get('s2', '')
    e2 = args.get('e2', '')
    
    fig = plt.figure()
    ax = plt.axes()
    
    
    xr = [0, length]
    yr = [1, 1]
    ax.plot(xr, yr);
    yr = [0, 0]
    ax.plot(xr, yr);
    yr = [-1, -1]
    ax.plot(xr, yr);

    x1 = []
    x2 = []
    x3 = []

    header_found = False
    num_words = 0
    for line in open(vcf_filename):
        words=line.split()
        if ((!header_found) && (words[1] == "POS")):
            header_found = True
            num_words = len(words)
        else:
            if (int(words[num_words-3]) > 0):
                x1.append(int(words[1]))
            if (int(words[num_words-2]) > 0):
                x2.append(int(words[1]))
            if (int(words[num_words-1]) > 0):
                x3.append(int(words[1]))

    y1 = [1 for tmp in x1]
    y2 = [0 for tmp in x2]
    y3 = [-1 for tmp in x3]
    ax.scatter(x1, y1, marker='o')
    ax.scatter(x2, y2, marker='o')
    ax.scatter(x3, y3, marker='o')

    ax.plot([s1,s1],[-10, 10],'k')
    ax.plot([s2,s2],[-10, 10],'k')
    ax.plot([e1,e1],[-10, 10],'k')
    ax.plot([e2,e2],[-10, 10],'k')

    plt.ylim([-10, 10])
    plt.show()
