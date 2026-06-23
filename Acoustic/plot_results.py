import os
import re
import numpy as np
import matplotlib.pyplot as plt

AU = 1.495978707e13
PC = 3.08567758e18


def parse(x):
    try:
        return float(x)
    except:
        return float(re.sub(r"([0-9])([+-])([0-9]+)$", r"\1E\2\3", x))


def load(fname):
    data = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                data.append([parse(t) for t in line.split()])
    return np.array(data)


def load_blocks(fname):
    blocks = []
    cur = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line == "":
                if cur:
                    blocks.append(np.array(cur))
                    cur = []
            elif not line.startswith("#"):
                cur.append([parse(t) for t in line.split()])
    if cur:
        blocks.append(np.array(cur))
    return blocks


# acoustic wave
def fig_acoustic():

    files = ["wave_t0.dat", "wave_thalf.dat", "wave_tT.dat", "wave_analytic.dat"]

    plt.figure(figsize=(10,4))

    for f in files:
        if os.path.isfile(f):
            d = load(f)
            x = d[:,0]/AU
            rho = d[:,1]
            rho0 = np.mean(rho)
            plt.plot(x,(rho-rho0)/rho0,label=f)

    plt.xlabel("x (AU)")
    plt.ylabel("drho/rho0")
    plt.legend()
    plt.title("acoustic wave")
    plt.savefig("fig1_acoustic.png")
    plt.show()


# jeans unstable
def fig_jeans_a():

    if os.path.isfile("jeans_growth_a.dat"):
        d = load("jeans_growth_a.dat")

        plt.figure()
        plt.semilogy(d[:,0],d[:,1],label="num")
        plt.semilogy(d[:,0],d[:,2],"--",label="ana")
        plt.xlabel("t")
        plt.ylabel("amp")
        plt.legend()
        plt.title("jeans unstable growth")
        plt.savefig("fig2_jeans_a.png")
        plt.show()


# jeans stable
def fig_jeans_b():

    if os.path.isfile("jeans_growth_b.dat"):
        d = load("jeans_growth_b.dat")

        A0 = d[0,2]

        plt.figure()
        plt.plot(d[:,0],d[:,1]/A0)
        plt.axhline(1,ls="--")
        plt.xlabel("t")
        plt.ylabel("amp/A")
        plt.title("jeans stable")
        plt.savefig("fig3_jeans_b.png")
        plt.show()



# marginal case
def fig_jeans_c():

    if not os.path.isfile("jeans_marginal.dat"):
        return

    d = load("jeans_marginal.dat")
    A0 = d[0,2]

    plt.figure()
    plt.plot(d[:,0],d[:,1]/A0)
    plt.axhline(1,ls="--")
    plt.xlabel("t")
    plt.ylabel("amp/A")
    plt.title("jeans marginal")
    plt.savefig("fig4_jeans_c.png")
    plt.show()



if __name__ == "__main__":

    fig_acoustic()
    fig_jeans_a()
    fig_jeans_b()
    fig_jeans_c()
    fig_dispersion()