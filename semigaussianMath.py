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

def findmGivennQSallenKey(n,Q,greaterThanOne=True):
    Q2 = Q**2
    firstTerm =  n/(2*Q2) - 1.
    secondTerm = sqrt(n*(n-4*Q2))/(2*Q2)
    if greaterThanOne:
      result = firstTerm+secondTerm
    else:
      result = firstTerm-secondTerm
    return result

def findmGivennQMFB(n,Q,greaterThanOne=True):
    Q2 = Q**2
    firstTerm =  (n-4*Q2)/(8*Q2) 
    secondTerm = sqrt(n**2-8*Q2*n)/(8*Q2)
    if greaterThanOne:
      result = firstTerm+secondTerm
    else:
      result = firstTerm-secondTerm
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
  e12Series = [1.0,1.2,1.5,1.8,2.2,2.7,3.3,3.9,4.7,5.6,6.8,8.2]
  e24Series = [1.0,1.1,1.2,1.3,1.5,1.6,1.8,2.0,2.2,2.4,2.7,3.0,3.3,3.6,3.9,4.3,4.7,5.1,5.6,6.2,6.8,7.5,8.2,9.1]
  print("{:1}  {:9}  {:9}  {:9}  {:9}  {:9}  {:9}".format("i","mSK","nSK","R2C2SK","mMFB","nMFB","R2C2MFB"))
  for i in range(1,len(w)):
      nSK = (Q[i-1]*(1+m))**2/m
      ax.semilogy(nSK,m,label="SK Term: {}".format(i))
      nMFB = (Q[i-1]*(1+2*m))**2/m
      ax.semilogy(nMFB,m,label="MFB Term: {}".format(i))
      nCalcSK = None
      mCalcSK = None
      #for nTry in e12Series:
      for nTry in e24Series:
        nCalcSK = nTry
        mCalcSK = findmGivennQSallenKey(nCalcSK,Q[i-1])
        if not isnan(mCalcSK):
          break
      ax.plot(nCalcSK,mCalcSK,label="SK n={:1.1f}, m={:.3f}".format(nCalcSK,mCalcSK),ls="",marker="o")
      R2C2SK = 1/(w[i]*sqrt(mCalcSK*nCalcSK))
      nCalcMFB = None
      mCalcMFB = None
      for nTry in e24Series:
        nCalcMFB = nTry
        mCalcMFB = findmGivennQMFB(nCalcMFB,Q[i-1])
        if not isnan(mCalcMFB):
          break
      R2C2MFB = 1/(w[i]*sqrt(mCalcMFB*nCalcMFB))
      ax.plot(nCalcMFB,mCalcMFB,label="MFB n={:1.1f}, m={:.3f}".format(nCalcMFB,mCalcMFB),ls="",marker="o")
      print("{:1}  {:9.7f}  {:9.7f}  {:9.7f}  {:9.7f}  {:9.7f}  {:9.7f}".format(i,mCalcSK,nCalcSK,R2C2SK,mCalcMFB,nCalcMFB,R2C2MFB))
  ax.legend(loc="best")
  fig.savefig("mVn{}.png".format(order))
  #print("{:1}  {:9}  {:9} ".format("i","m1","m2"))
  #for i in range(1,len(w)):
  #    m1,m2 = mFornAndQ(10.,Q[i-1])
  #    print("{:1}  {:9.4f}  {:9.4f} ".format(i,m1,m2))

  

