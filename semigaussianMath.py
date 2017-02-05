#!/usr/bin/env python3

from numpy import *
from matplotlib import pyplot as plt

def lowPassH(s,w):
    """
    s = laplace variable
    w = angular frequency of cutoff
    """
    return 1.0/(s+w)

def quadraticLowPassH(s,Q,w):
    """
    s = laplace variable
    Q is Q-factor
    w = angular frequency of cutoff
    """
    w2 = w**2
    return w2/(s**2+s*w/Q + w2)

def semiGaussianLowPassH(s,wList,QList):
    if (len(wList) + len(QList) % 2) == 0:
      raise NotImplementedError("Even order not implemented")
    else:
      result = lowPassH(s,wList[0])
      for i in range(1,len(wList)):
        result *= quadraticLowPassH(s,wList[i],QList[i-1])
    return result

# A
real5 = [1.4766878,1.4166647,1.2036832]
real7 = [1.6610245,1.6229725,1.4949993,1.2344141]
# W
imag5 = [0.5978596,1.2994843]
imag7 = [0.5007975,1.0454546,1.7113028]

print("w and alpha are per sigma")
for real, imag in [(real5,imag5),(real7,imag7)]:
  A = array(real)
  W = array(imag)

  # omega^2 = (A^2 + W^2) / sigma^2
  w = sqrt(A[:-1]**2 + W**2)
  
  # alpha = A / sigma
  alpha = A[:-1]
  
  # Q = omega / alpha / 2
  Q = w/alpha/2

  # put A0 on front of w
  w = insert(w,0,A[0])

  order = len(real)+len(imag)
  print("Order: ",order)
  print("{:1}  {:9}  {:9}".format("i","w","Q"))
  #print("{:1}  {:9.7f}".format(0,A[0]))
  for i in range(len(w)):
      if i > 0:
        print("{:1}  {:9.7f}  {:9.7f}".format(i,w[i],Q[i-1]))
      else:
        print("{:1}  {:9.7f}  ".format(i,w[i]))

  fig, ax = plt.subplots()
  ax.set_xlabel("n=C1/C2")
  ax.set_ylabel("m=R1/R2")
  ax.set_title("Allowed Values of m given Q and n")
  m = logspace(-1,1)
  for i in range(1,len(w)):
      nCalc = (Q[i-1]*(1+m))**2/m
      ax.semilogy(nCalc,m,label="Term: {}".format(i))
  ax.legend(loc="best")
  fig.savefig("mVn{}.png".format(order))
  #print("{:1}  {:9}  {:9} ".format("i","m1","m2"))
  #for i in range(1,len(w)):
  #    m1,m2 = mFornAndQ(10.,Q[i-1])
  #    print("{:1}  {:9.4f}  {:9.4f} ".format(i,m1,m2))

  

