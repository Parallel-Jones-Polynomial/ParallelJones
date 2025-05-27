import math
import numpy as np
from copy import deepcopy
import time
from itertools import permutations

def generate_tuples(N):
	M=set([])
	for X in range(N+1):
		indices = [1] * X + [-1] * (N - X)
		tuples = set(permutations(indices, N))
		M=M.union(tuples)
	return M

def print_binary_combinations(n):
	# Loop through all numbers from 0 to 2^n - 1
	L=[]
	for i in range(1<<n):
   	# Convert the current number to a binary string of length n
		binary_str = format(i, '0' + str(n) + 'b')
		L.append(binary_str)
	return L	

import itertools
def prod_iter(N):
	L=[]
	W = itertools.product([-1,1], repeat=N)
	for w in W:
		L.append(w)
	return L	


# Picks a random point from the surface of a unit sphere.
def get_random_proj():
	b = 0
	while np.sum(b*b)>1 or np.sum(b*b)<1e-4:
		b = 2*np.random.uniform(size=3)-1
	return b/np.sqrt(np.sum(b*b))


# Chooses n = samples points from the surface of a unit sphere.
def fibonacci_sphere(samples):
	points = []
	phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
	for i in range(samples):
		y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
		radius = math.sqrt(1 - y * y)  # radius at y
		theta = phi * i  # golden angle increment
		x = math.cos(theta) * radius
		z = math.sin(theta) * radius
		points.append((x, y, z))
	return points


def get_two_vec(proj_vec): 
   a = np.zeros([3])
   proj_vec = np.array(proj_vec)
   while np.sum(a*a) < 1e-7:
     a = np.random.normal(loc=-10.0, scale=1e-7, size=[3])
     a = a - proj_vec*np.sum(a*proj_vec)/np.sum(proj_vec*proj_vec)
   a /= np.sqrt(np.sum(a*a))
   b = np.array([
          a[1]*proj_vec[2]-a[2]*proj_vec[1], 
           -a[0]*proj_vec[2]+a[2]*proj_vec[0],
            a[0]*proj_vec[1]-a[1]*proj_vec[0]
            ])   #cross product with a and proj
               #orthogonality of a and b can be checked using np.dot(a,b)
   #print('a, b, dot cross',np.array(a),np.array(b),np.dot(proj_vec, np.cross(a,b)))
   if np.dot(proj_vec, np.cross(a,b))>0:
      return a, b
   elif np.dot(proj_vec, np.cross(a,b))<0:
      return b, a

## General method to find connected components of a graph
def getRoots(aNeigh):
   def findRoot(aNode,aRoot):
      while aNode != aRoot[aNode][0]:
        aNode = aRoot[aNode][0]
      return (aNode,aRoot[aNode][1])
   myRoot = {} 
   for myNode in aNeigh.keys():
      myRoot[myNode] = (myNode,0)  
   for myI in aNeigh: 
      for myJ in aNeigh[myI]: 
        (myRoot_myI,myDepthMyI) = findRoot(myI,myRoot) 
        (myRoot_myJ,myDepthMyJ) = findRoot(myJ,myRoot) 
        if myRoot_myI != myRoot_myJ: 
           myMin = myRoot_myI
           myMax = myRoot_myJ 
           if  myDepthMyI > myDepthMyJ: 
             myMin = myRoot_myJ
             myMax = myRoot_myI
           myRoot[myMax] = (myMax,max(myRoot[myMin][1]+1,myRoot[myMax][1]))
           myRoot[myMin] = (myRoot[myMax][0],-1) 
   myToRet = {}
   for myI in aNeigh: 
      if myRoot[myI][0] == myI:
        myToRet[myI] = []
   for myI in aNeigh: 
      myToRet[findRoot(myI,myRoot)[0]].append(myI) 
   return myToRet  

def Loops(inds,clos_perm):
	#print('Loopy',inds)
	#print('clos_perm inside Loops',clos_perm)
	## Take into account relation between indices due to closure permutation.
	working_inds=deepcopy(inds)
	try:
		for i in clos_perm:
			#print('iii',i)
			try:
				working_inds[i].append(clos_perm[i])
			except:
				pass
	except:
		pass
	S={} #preimages
	for i in working_inds:
		for j in working_inds[i]:
			try:
				S[j]+=[i]
			except:
				S[j]=[i]		
	## Neighbors of any i (a dictionary)
	NGB=deepcopy(working_inds)
	for i in S:
		if i in NGB:
			NGB[i]+=S[i]
		else:
			NGB[i]=S[i]
	## Partition vertices into mutually disjoint sets.
	L_raw=getRoots(NGB)
	L=[set(L_raw[i]) for i in L_raw]
	## Ensure that elts in L are mutually disjoint
	#print('NGB',NGB)
	#print('L',L)			
	return NGB, L

''' CRAMER'S RULE 2x2'''
def det2x2(A):
	assert A.shape == (2,2)
	return A[0][0]*A[1][1] - A[0][1]*A[1][0]

def Cramer(A):
	assert A.shape == (2,3)
	D = det2x2(A[:,:2])
	if D == 0:
		return
	Dx = det2x2(A[:,[2,1]])
	Dy = det2x2(A[:,[0,2]])
	return Dx*1.0/D, Dy*1.0/D



# Reidemeister MOVES
def RM1(bool_mask,inds):
	#print('total number of crossings before RM 1 :',np.count_nonzero(bool_mask)/2.)
	Bin=[]
	for a in inds:
		if a not in Bin:
			if np.any(bool_mask[a, :]):
				flag=0
				b=inds[a]
				#print(a,b)
				while flag==0:
					if np.any(bool_mask[b, :]):
						if bool_mask[b,a]==True:
							bool_mask[a,b]=False
							bool_mask[b,a]=False
							Bin.append(a)
							Bin.append(b)
							flag=1
						else:
							Bin.append(a)	
							flag=1
					else:
						try:
							b=inds[b]
						except:
							flag=1
	#print('total number of crossings after RM 1 :',np.count_nonzero(bool_mask)/2.)
	return bool_mask

def RM2(bool_mask, over_or_under, inds):
   Bin = []
   #print('total number of crossings before RM 2 :',np.count_nonzero(bool_mask)/2.)
   #print(bool_mask, '\n', over_or_under, '\n', inds)
   for a in inds:
      if a not in Bin:
        b = inds[a]
        if np.any(bool_mask[a, :]):
           c = np.argmax(bool_mask[a, :])
           d = inds[c]
           try:
             if bool_mask[b, d] == True:
               print('elementary RM2 happens')
               if over_or_under[a, c] == over_or_under[b, d]:
                   bool_mask[a, c] = False
                   bool_mask[c, a] = False
                   bool_mask[b, d] = False
                   bool_mask[d, b] = False
                   Bin.append(a)
                   Bin.append(b)
                   Bin.append(c)
                   Bin.append(d)
           except:
                flag=0
                while flag==0:
                    if np.any(bool_mask[b, :]):
                        g=np.argmax(bool_mask[b, :])
                        h=inds[g]
                        if find_conn(d, g, inds, bool_mask)==1:
                            if g not in Bin and  over_or_under[a, c] == over_or_under[b, g]:
                                bool_mask[a, c] = False
                                bool_mask[c, a] = False
                                bool_mask[b, g] = False
                                bool_mask[g, d] = False
                                Bin.extend([a, c, b, g])  
                        if find_conn(h, c, inds, bool_mask)==1:
                            if h not in Bin and  over_or_under[g, b] == over_or_under[c, a]:
                                bool_mask[a, c] = False
                                bool_mask[c, a] = False
                                bool_mask[b, g] = False
                                bool_mask[g, d] = False
                                Bin.extend([a, c, b, g])    
                    else:
                        try:
                            b=inds[b]
                        except:
                            flag=1
		   

   #print('total number of crossings after RM 2 :',np.count_nonzero(bool_mask)/2.)
   return bool_mask

def find_conn(d, g, inds, bool_mask):
    V = 0
    flag = 0
    while flag == 0:
        try:
            if d == g:
                V = 1
                flag = 1
            else:
                try:
                    d = inds[d]
                except:
                    flag = 1
        except:
            flag = 1
    return V

def simplification(BM,over_or_under,inds):
	print('total number of crossings before RM :',np.count_nonzero(BM)/2.)
	#print('yep')
	#BM=RM3(BM, over_or_under,inds)
	#BM=RM2(BM, over_or_under,inds)
	#print(BM)
	
	#BM=RM1(BM,inds)
	#BM=RM2(BM,over_or_under,inds)
	# for i in range(5):
	# 	print(i,'-----')
	# 	BM=RM1(BM,inds)	
	for i in range(1):
		print(i,'******')
		BM=RM1(BM,inds)	
		BM=RM2(BM,over_or_under,inds)
		BM=RM1(BM,inds)	
	#for i in range(5):
	#	BM=RM1(BM,inds)	
	#print(BM)
	print('total number of crossings after RM :',np.count_nonzero(BM)/2.)
	return BM		

def max_len(Ch):
	l=0
	for ch in Ch:
		l=max([l,len(ch)])
	return l

def J_mult(P,Q):
	Z={}
	for i in P:
		for j in Q: 
			try:
				Z[i+j]+=P[i]*Q[j]
			except:
				Z[i+j]=P[i]*Q[j]
	return Z

def J_add(P,Q):
	Z=P.copy()
	for j in Q: 
		try:
			Z[j]+=Q[j]
		except:
			Z[j]=Q[j]
	return Z

def J_smult(a,P):
	Z={}
	for i in P:
		Z[i]=a*P[i]
	return Z
	
def dfactor(N):
	dpoly={0:1}
	for i in range(N):
		dpoly=J_mult(dpoly,{-2:-1,2:-1})
	return dpoly


