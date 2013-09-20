import os

def main():
    dir = os.getcwd()
    if dir[-5:] != "Tests":
        os.chdir("Tests")
    names = os.listdir(os.curdir)
    for name in names:
        if name[:5] == "test_" and name[-3:] == ".py":
            print name
            os.system("python " + name)

if __name__ == "__main__":
    main()
