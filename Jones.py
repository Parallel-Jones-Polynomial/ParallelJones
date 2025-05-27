'''
JONES POLYNOMIAL PROGRAM
=========================
Given the coordinates of a collection of open/closed curves in 3-space
the following program calculates its Jones polynomial.
'''
import math
import numpy as np
from copy import deepcopy
import multiprocessing
import time
from functions import *
from sys import getsizeof 


''' 
get_jones_poly : 
coords signifies the co-ordinates of points on the curve in 3space
proj_vector signifies a particular direction of projection
inds signifies the connections among the points in the curve in concern.
'''
"""Calculate the Jones polynomial wrt a projection vector."""
def get_jones_poly(coords, proj_vector, Inds):
    #print('Inds inside get jones poly:',Inds)
    inds=Inds[0]
    clos_perm=Inds[1]
    #print('inds',inds)
    #print('clos_perm',clos_perm)
    #print("#2")
    ## Generate two vectors from the projection vector
    x, y = get_two_vec(proj_vector)
    ## Project the coordinates onto the plane spanned by the two vectors
    W=np.array([x,y,proj_vector])
    proj_new={}
    depth={}
    for i in coords:
        v=coords[i]
        v_new=np.linalg.solve(W, v) 
        proj_new[i]=[v_new[0], v_new[1]]
        depth[i]=v_new[2]
    ## Uncomment to view (x,y) coordinates along proj vec   
    #print('New',proj_new)
    ## Uncomment to view depth aka new z coordinate (along proj vec)
    #print('depth',depth)    
    ## Check if the projection is irregular i.e. if 2 or more vertices are overlapping.
    Check_proj = np.unique([str(proj_new[i]) for i in proj_new], return_counts=True)
    if max(Check_proj[1]) >= 2:
        #print('irreg points',max(Check_proj[1]))
        return None
    ## Generate crossing matrices and associated properties.
    bool_mask, over_or_under, right_or_left, inds = crossing_matrices_bool(proj_new, inds, depth)
    #print('inds after addition', inds)
    ## Calculate Writhe.
    Wr = get_writhe(bool_mask, over_or_under, right_or_left)
    ## Uncomment to print the writhe and no. of crossings wrt the given projection.
    #print('Writhe before RM', Wr)
    #print('total number of crossings before RM :',np.count_nonzero(bool_mask)/2.)
    ## Do Reidemeister Moves
    if np.count_nonzero(bool_mask)/2.>=0.:
        bool_mask=simplification(bool_mask, over_or_under,inds)
        Wr = get_writhe(bool_mask, over_or_under, right_or_left)
        ## Uncomment to print the writhe and no. of crossings wrt the given projection.
        #print('Writhe after RM', Wr)
        #print('total number of crossings after RM :',np.count_nonzero(bool_mask)/2.)
        ## Identify very high crossing configurations
        if np.count_nonzero(bool_mask)/2.>36:
            return None
            #return np.count_nonzero(bool_mask)/2.
        else:
            pass
    ## Parallel Computation : Via Splitting and Gluing
    if parallel==1:
        print('More',np.count_nonzero(bool_mask)/2. )
        #print('enters splitting')
        ## SPLITTING OF THE KNOT INTO 2 LINKOIDS L1 and L2
        ## Arc defining L1
        Arc=find_arc(bool_mask,inds)
        #print('Arcs I and J',Arc)
        L1_arcs=Arc[0]+Arc[1]
        inds_L1 = {k: v for k, v in inds.items() if k in L1_arcs}
        inds_L2 = {k: v for k, v in inds.items() if k not in L1_arcs}
        #print('inds_L1',inds_L1)
        #print('inds_L2',inds_L2)
        L2_arcs = list(inds_L2.keys())
        L2_arcs.sort()
        #print('L1_arcs',L1_arcs)
        #print('L2_arcs',L2_arcs)
        BM_L1= np.copy(bool_mask)
        OU_L1= np.copy(over_or_under)
        RL_L1= np.copy(right_or_left)
        for i in L2_arcs:
            BM_L1[i,:]=False
            BM_L1[:,i]=False
            OU_L1[i,:]=False
            OU_L1[:,i]=False
            RL_L1[i,:]=False
            RL_L1[:,i]=False
        BM_L2= np.copy(bool_mask)
        OU_L2= np.copy(over_or_under)
        RL_L2= np.copy(right_or_left)
        for i in L1_arcs:
            BM_L2[i,:]=False
            BM_L2[:,i]=False
            OU_L2[i,:]=False
            OU_L2[:,i]=False
            RL_L2[i,:]=False
            RL_L2[:,i]=False
        ## First bring inds_L1 and inds_L2 to desired format.
        modf_inds_L1 = {}
        modf_inds_L2 = {}
        for i in inds_L1:
            modf_inds_L1[i]=[inds_L1[i]]
        for i in inds_L2:
            modf_inds_L2[i]=[inds_L2[i]]
        states1=[]
        states2=[]
        arg_1=BM_L1, OU_L1, RL_L1, inds_L1,{}, modf_inds_L1,states1,list(np.copy(states1))
        arg_2=BM_L2, OU_L2, RL_L2, inds_L2,{}, modf_inds_L2,states2,list(np.copy(states2))
        #print('A', len(arg_1), len(arg_2))
        ## Return the Jones polynomial in terms of the bracket polynomial states of the L1 and L2 linkoids.
        SJ=split_jones(arg_1,arg_2,clos_perm,Wr)
        #print('Bkt with splitting',SJ)
        return SJ
    ## Serial Computation
    else:
        print('Less', np.count_nonzero(bool_mask)/2. )
        ## Compute the Bkt and Jones without splitting.
        states=[]
        modf_inds = {}
        for i in inds:
            modf_inds[i]=[inds[i]]
        #print('Inds before get PP',Inds)
        start_time=time.time()
        L_bkt = get_partial_poly(bool_mask, over_or_under, right_or_left,inds,clos_perm, modf_inds,states,list(np.copy(states)))
        ## Uncomment to print the bracket polynomials of the L1 and L2 linkoids.
        print("L_bkt size", getsizeof(L_bkt))
        Bkt_poly={0:0}
        for K in L_bkt:
            Bkt_poly=J_add(Bkt_poly,J_mult({K[0]:1},dfactor(len(K[1])-1)))
        end_time=time.time()
        print('Time taken serial :', end_time-start_time)
        #print('Bkt w/o splitting',Bkt_poly)
        #print('Jones w/o splitting',J_mult(Bkt_poly,{3*(-Wr):(-1)**abs(Wr)}))
        return J_mult(Bkt_poly,{3*(-Wr):(-1)**abs(Wr)})
        
"""Create crossing matrices, new vertices and other properties """
def crossing_matrices_bool(proj, inds, depth):
    #print("#3")
    ## Create ordered edges from inds dict
    edges={}
    # ct=0
    # for i in inds.keys():
    #     edges[ct]=(i,inds[i])
    #     ct+=1 
    ct=0
    for i in inds.keys():
        edges[i]=(i,inds[i])
        ct+=1 
    ## Uncomment to view original inds, proj and edges.  
    #print('inds1',inds)    
    #print('edges1',edges)
    #print('proj1',proj)   
    #print('depths1',depth) 
    ## Create 4 types of matrices :
    ## bool_mask : a_ij = True (False) if edge i and edge j have (does not have) a crossing. 
    ## over_or_under : a_ij = True (False) if edge i lies above (does not lie above) edge j. 
    ## right_or_left : a_ij = True (False) if cross product: edge i x edge j is positive (negative). (Right hand thumb rule).
    ## u :  a_ij = coords of intersection point (if any) b/w edge i and edge j.
    modf_inds0 = {}
    for i in inds:
        modf_inds0[i]=[inds[i]]
    #print('modf_inds0 1 :', modf_inds0)
    Loops_as_sets=Loops(modf_inds0,{})[1]
    #print('Loops_as_sets 1 :', Loops(modf_inds0,{})[1])
    total_inds=[]
    for elt in Loops_as_sets:
        total_inds+=list(elt)
    #print('total inds 1 :',total_inds)
    I_inds=len(total_inds)
    bool_mask = np.zeros((I_inds,I_inds), dtype=bool)
    over_or_under = np.zeros((I_inds,I_inds), dtype=bool)
    right_or_left = np.zeros((I_inds,I_inds), dtype=bool)
    u = np.zeros((I_inds,I_inds), dtype=object)
    #print('BM first', bool_mask)
    ## dif : stores vector rep for any edge, (p1,p2), in terms of coords of the endpoints : vec(p2)-vec(p1).
    ## flip neg : similar to dif but replace [x,y] -> [y,-x] (flipped and negative). Used in Cramer rule later.
    dif={}
    flip_neg={}
    for e in edges.keys():
        #print('e',e,'edges[e]',edges[e])
        #print(proj[edges[e][0]],proj[edges[e][1]])
        dif[e]=[proj[edges[e][1]][0]-proj[edges[e][0]][0],proj[edges[e][1]][1]-proj[edges[e][0]][1]]
        flip_neg[e]=[proj[edges[e][1]][1]-proj[edges[e][0]][1],-proj[edges[e][1]][0]+proj[edges[e][0]][0]]
    ## Uncomment to view dif and flip_neg dictionaries.    
    #print('diff', dif)
    #print('flip neg diff', flip_neg)
    ## Iterate through each pair of edges to find intersection point(s) and update the matrices.
    ## intersection pt vs edges which contain it.
    intp_vs_edges = {} 
    num_edges=edges.keys()
    for ii in range(len(num_edges)):
        for jj in range( ii+1,len(num_edges)):
            i=list(num_edges)[ii]
            j=list(num_edges)[jj]
            #print('edgesss',ii,jj,i,j)
            #print('Diffs',dif[i],dif[j])
            ## cross product +ve or -ve
            rl_ij=np.cross(dif[i],dif[j])
            ## Compute terms for the linear system to find intersections.
            pi=proj[edges[i][0]][0]*proj[edges[i][1]][1]-proj[edges[i][0]][1]*proj[edges[i][1]][0]
            pj=proj[edges[j][0]][0]*proj[edges[j][1]][1]-proj[edges[j][0]][1]*proj[edges[j][1]][0]
            Aij=np.array([[flip_neg[i][0],flip_neg[i][1],pi],[flip_neg[j][0],flip_neg[j][1],pj]])
            ## Solve the linear system to find intersections using Cramer's rule.
            result=Cramer(Aij)
            try:
                int_pt=[result[0], result[1]]
            except:
                int_pt=None    
            ## Uncomment to view intersection point between i and j.    
            #print('intersection point',i, j, int_pt)
            ## If an intersection is found, matrices are updated as follows: 
            if int_pt != None:
                #print('check Right or Left:',i,j,'R or L',rl_ij)
                if rl_ij==0:
                    return None
                ## Parameterize intersection points in terms of s and t.
                xdif_i=proj[edges[i][1]][0]-proj[edges[i][0]][0]
                xdif_j=proj[edges[j][1]][0]-proj[edges[j][0]][0]
                if xdif_i!=0:
                    t=(int_pt[0]-proj[edges[i][0]][0])/(proj[edges[i][1]][0]-proj[edges[i][0]][0])
                else:
                    t=(int_pt[1]-proj[edges[i][0]][1])/(proj[edges[i][1]][1]-proj[edges[i][0]][1])
                if xdif_j!=0:
                    s=(int_pt[0]-proj[edges[j][0]][0])/(proj[edges[j][1]][0]-proj[edges[j][0]][0])
                else:
                    s=(int_pt[1]-proj[edges[j][0]][1])/(proj[edges[j][1]][1]-proj[edges[j][0]][1])
                ## Set a numerical threshold for zero
                zero_val=1e-4
                ## Uncomment to view intersection point between i and j.    
                #print(i,j,'intersection point', int_pt, t, s)
                ## Check if intersection point is a crossing i.e. lies within the edge segments.
                if zero_val<t<1-zero_val and zero_val<s<1-zero_val: 
                    #print(i,j,' good intersection point', int_pt, t, s)
                    bool_mask[i,j]=True
                    bool_mask[j,i]=True
                    intp_vs_edges[str(int_pt)]=[i,j]
                    ## Record which edge is over the other at the intersection point in over_or_under matrix.
                    zi=depth[edges[i][0]]+t*(depth[edges[i][1]]-depth[edges[i][0]])
                    zj=depth[edges[j][0]]+s*(depth[edges[j][1]]-depth[edges[j][0]])
                    #print('edge1:',i, 'edge2:', j, 'R or L:',rl_ij,'O or U:',zi-zj)
                    if zi>zj:
                        ## edge i lies over edge j.
                        over_or_under[i,j]=True 
                        #over_or_under[j,i]=True 
                    else:
                        ## edge j lies over edge i
                        over_or_under[j,i]=True 
                        #over_or_under[i,j]=True
                    if rl_ij>0:
                        right_or_left[i,j]=True
                        #right_or_left[j,i]=False
                    elif rl_ij<0:
                        right_or_left[j,i]=True 
                        #right_or_left[i,j]=False
                    ## Record the intersection point.
                    u[i,j]=int_pt
                    u[j,i]=int_pt  
    ## Uncomment to view matrices before new vertices are introduced. 
    #print('original indices',inds)                
    #print('bool_mask before new vertices:', bool_mask[25])
    #print('over_or_under before new vertices:',over_or_under)
    #print('right_or_left before new vertices:',right_or_left)
    #print('intersection points matrix before new vertices:', u)

    ## Introduction of new vertices by creating subdivision of edges wherever one edge 
    ## participates in > 1 crossings. Runs only once.
    ## intersections_edge stores edge vs. crossings on the edge.
    intersections_edge={}
    for i in range(len(total_inds)):
        if np.count_nonzero(bool_mask[i,:])>=1:
            pts_on_edge=[]
            for j in range(len(total_inds)):
                if bool_mask[i,j]==True:
                    pts_on_edge.append(u[i,j])                  
            intersections_edge[i]=pts_on_edge
    ## Uncomment to view intersections_edge dictionary.        
    #print('New points', intersections_edge) 
    #for i in intersections_edge:
        #print(i, intersections_edge[i])
    #print('-------')
    ## v denotes new index label.        
    v=len(proj.keys()) 
    ## intersection pt vs new edges which contain it.
    intp_vs_newedges = {}       
    ## i denotes an orginial edge that is to be subdivided.
    for i in intersections_edge:
        e=edges[i]
        int_points=intersections_edge[i]
        ## Create a dictionary to store absolute magnitude of 
        ## displacement from the starting pt of the edge vs. the intersection pt in concern. 
        mag_vs_pt = {}
        for elt in int_points:
            mag_vs_pt[math.dist(elt, proj[e[0]])] = elt
        ## Uncomment to view mag_vs_pt dictionary.    
        #print('Mag',e, mag_vs_pt)    
        ## Arrange the intersection points in ascending order in terms of their displacement magnitudes.
        mag_ord_keys = list(mag_vs_pt.keys())
        mag_ord_keys.sort()
        ## ini_v is the label of the stating point of an edge.
        ini_v=e[0]
        for j in range(len(mag_ord_keys)):
            ## For every pair of consecutive intersection points along the edge,
            ## create a new point at their midpoint. and add it to the indslist of vertices.
            p1 = mag_vs_pt[mag_ord_keys[j]]
            try:
                p2 =  mag_vs_pt[mag_ord_keys[j+1]]
                new_p = [(p1[0] + p2[0]) / 2., (p1[1] + p2[1]) / 2.]
                ## add new point to proj, inds and edges dictionaries. 
                try:
                    intp_vs_newedges[str(p1)].append(ini_v)
                except:
                    intp_vs_newedges[str(p1)]=[ini_v]    
                proj[v]=new_p
                edges[ini_v]=(ini_v,v)
                inds[ini_v]=v
                ini_v=v
                v+=1
            except:
                p2 = proj[e[1]]    
                try:
                    intp_vs_newedges[str(p1)].append(ini_v)
                except:
                    intp_vs_newedges[str(p1)]=[ini_v]
                edges[ini_v]=(ini_v,e[1])
                inds[ini_v]=e[1]
    ## Uncomment to print intersection pt vs participating edges (old and new)            
    #print('IVE', intp_vs_edges)            
    #print('IVNE',intp_vs_newedges)
    ## Uncomment to view the updated indices, proj and edges.
    #print('inds2', inds)
    #print('edges2',edges)
    #print('proj',proj)
    ## Update the matrices to account for the newly introduced vertices.
    ## BM_2, OU_2 and RL_2 are expanded to make space for subdivided edges.
    modf_inds1 = {}
    for i in inds:
        modf_inds1[i]=[inds[i]]
    Loops_as_sets1=Loops(modf_inds1,{})[1]
    total_inds1=[]
    for elt in Loops_as_sets1:
        total_inds1+=list(elt)
    J_inds=len(total_inds1)
    BM_2= np.zeros((J_inds,J_inds), dtype=bool)
    OU_2= np.zeros((J_inds,J_inds), dtype=bool)
    RL_2= np.zeros((J_inds,J_inds), dtype=bool)
    for pt in intp_vs_newedges:
        #print('PTTT',pt,intp_vs_newedges[pt])
        i=intp_vs_newedges[pt][0]
        j=intp_vs_newedges[pt][1]
        BM_2[i,j]=True
        BM_2[j,i]=True
        i0=intp_vs_edges[pt][0]
        j0=intp_vs_edges[pt][1]
        #print('ptsss',i,j,i0,j0)
        OU_2[i,j]=np.copy(over_or_under)[i0,j0]
        OU_2[j,i]=np.copy(over_or_under)[j0,i0]
        RL_2[i,j]=np.copy(right_or_left)[i0,j0]
        RL_2[j,i]=np.copy(right_or_left)[j0,i0] 
    #print('BM_2',BM_2[25])
    #print('new_inds',inds)
    return BM_2, OU_2, RL_2, inds

"""Calculate BRACKET POLYNOMIAL using crossing matrices and associated properties"""
def get_partial_poly(bool_mask, over_or_under, right_or_left, inds, clos_perm, modf_inds,states,S):
    #print("#4")
    #print('States',states)
    ## Check any entry in bool_mask is True
    #print("L_bkt size", getsizeof(L_bkt))
    if np.any(bool_mask):
        #print('SS',states)
        edge1 = np.argmax(np.any(bool_mask, 0))
        edge2 = np.argmax(bool_mask[edge1,:])
        #print('indices', modf_inds)
        #print('edge1, edge2 :', edge1, edge2) 

        bool_mask1 = np.copy(bool_mask)
        bool_mask1[edge1, edge2] = False
        bool_mask1[edge2, edge1] = False
        bool_mask2 = np.copy(bool_mask1) 
        modf_inds[edge1].remove(inds[edge1])
        modf_inds[edge2].remove(inds[edge2])

        modf_inds1=deepcopy(modf_inds)
        modf_inds2=deepcopy(modf_inds)
    
        states1=list(np.copy(states))
        states2=list(np.copy(states))
        #print('edges',edge1,edge2,'i over j', over_or_under[edge1, edge2], 'i right j', right_or_left[edge1, edge2])
        if over_or_under[edge1, edge2] == right_or_left[edge1, edge2]:
            #print('pos',edge1,edge2)
            modf_inds1[edge1].append(inds[edge2])
            modf_inds1[edge2].append(inds[edge1])
            modf_inds2[edge1].append(edge2)
            try:
                modf_inds2[inds[edge1]].append(inds[edge2]) 
            except:
                modf_inds2[inds[edge1]]=[inds[edge2]]
            #print('M1,M2', modf_inds1, modf_inds2)
            #states1.append(1)
            #states2.append(-1)
        else:
            #print('neg',edge1,edge2)
            modf_inds1[edge1].append(edge2) 
            try:
                modf_inds1[inds[edge1]].append(inds[edge2])
            except:
                modf_inds1[inds[edge1]]=[inds[edge2]]
            
            modf_inds2[edge1].append(inds[edge2])
            modf_inds2[edge2].append(inds[edge1])
            #print('M1,M2', modf_inds1, modf_inds2)
            #states1.append(-1)
            #states2.append(1)
        states1.append(1)
        states2.append(-1)
        #print('States', states1)
        #print('pp1')
        get_partial_poly(bool_mask1, np.copy(over_or_under), np.copy(right_or_left), inds, clos_perm, modf_inds1,states1,S)
        #print('pp2')
        get_partial_poly(bool_mask2, np.copy(over_or_under), np.copy(right_or_left), inds, clos_perm, modf_inds2,states2,S)
    else:
        #print('clos perm before USE', clos_perm)
        S.append((sum(states),Loops(modf_inds,clos_perm)[1]))
    #print('SSSS',S)
    return S

''' WRITHE '''
def get_writhe(bool_mask, over_or_under, right_or_left):
    #print("#5")
    BM=np.copy(bool_mask)
    #print('pos cross',np.int32(over_or_under == right_or_left)[BM])
    #print('neg cross',np.int32(over_or_under != right_or_left)[BM])
    Wr=(np.sum(np.int32(over_or_under == right_or_left)[BM])-np.sum(np.int32(over_or_under != right_or_left)[BM]))/2.
    #print('Writhe Value',Wr)
    return Wr

## Find the arc that forms the backbone of L1.
## I is the backbone arc and J stores the auxillary arcs.
def find_arc(bool_mask,inds):
    #print('Inds fA',inds)
    C=int((np.count_nonzero(bool_mask)/2.)/2.)
    if C==0:
        C=1
    #print('cross',C)
    c=0
    i=0
    I=[]
    J=[]
    #print(' MY inds',inds)
    while c<C:
        I.append(i)
        #print('key check',i,inds)
        #c+=1
        j=np.argmax(bool_mask[i,:])
        if bool_mask[i,j]==True:
            #print('found a crossing')
            J.append(j)
            c+=1
        #print('inds Arc',inds)
        #print('I and J', I, J)
        try:
            i=inds[i]
        except:
            i = i+1
    #print('IZZZZ',I,J)
    return I,J

def form_indices(Ch,closed):
    ini_inds0={}
    ini = 0
    ## Index labels vs. coordinates in 3-space
    coord_dict = {}
    clos_perm={}
    #print('clos1',clos_perm)
    ## Generate indices for each component in Ch
    for ch in Ch:
        ## Uncomment to print each component.
        #print('component',ch)
        for i in range(len(ch)-1):
            ini_inds0[i+ini]=ini+i+1
            coord_dict[i+ini]=list(ch[i])
        #ini_inds[ini+len(ch)-1]=ini 
        #coord_dict[ini+len(ch)-1]=list(ch[len(ch)-1])
        if closed==1:
           ini_inds0[ini+len(ch)-1]=ini 
           coord_dict[ini+len(ch)-1]=list(ch[len(ch)-1])
        else:
           clos_perm[ini+len(ch)-1]=ini
           clos_perm[ini]=ini+len(ch)-1
           coord_dict[ini+len(ch)-1]=list(ch[len(ch)-1])    
        ini+=len(ch)   
    #print('clos2',clos_perm)
    ini_inds=[ini_inds0,clos_perm]
    #print('ini_inds',ini_inds)
    return coord_dict,ini_inds

def process_projection(knot_data):
    # Jones polynomial along a specific projection vector
    #print('KD',knot_data)
    rj = get_jones_poly(*knot_data)
    #print('rj',rj)
    return rj


def two_knot_jones(K,n):
    if __name__ == '__main__':
        points = []
        if n > 1:
            points = fibonacci_sphere(n)
        else:  
            points = [np.array([0., 0., 1.])]
        inputs=[(form_indices(K[i],closed)[0], p, form_indices(K[i],closed)[1]) for i in range(len(K)) for p in points]
        #print('ip',inputs)

        with multiprocessing.Pool() as pool:
            results = pool.map(process_projection, inputs)
            #results1 = pool.map(process_projection, inputs1)
            #results2 = pool.map(process_projection, inputs2)
            #results1.join()
        JPOLY = results
        #print('len JPOLY',len(JPOLY))
        #print('JPOLY',JPOLY)


def process_partial_poly(knot_data):
    #print('entered process partial poly')
    rj = get_partial_poly(*knot_data)
    #print('rj process partial poly:', rj)
    return rj

from multiprocessing import Pool

def evaluate_with_f(buf):
    # Create tuples with 2 elements from the buffer
    if len(buf) == 2:
        return tuple(buf)
    else:
        raise ValueError("error.")

def process_chunk(args):
    W1_chunk, W2 = args
    buf = []
    R = []
    #print('process')
    for i in range(len(W1_chunk)):
        #print(W1_chunk[i])
        # Copy element i of W1_chunk to the head of buf
        buf.append(W1_chunk[i])
        for k in range(len(W2)):
            # Append the kth element of W2 to buf
            buf.append(W2[k])

            # Evaluate one element of W with f(buf)
            result = evaluate_with_f(buf)
            #print('res',result)
            R.append(result)
            #print(f"Result: {result}")

            # Reset buf for the next iteration
            buf.pop()

        # Reset buf for the next iteration
        buf = []
    #print('R',R)
    return R

def mult_2_elements(L_M):
    L, M, clos_perm = L_M
    #clos_perm=L_M[1]
    #print('PPPPPP',L,M)
    Af = L[0] + M[0]
    S = {}
    for i in range(len(L[1])):
        S[i] = L[1][i]
    for j in range(len(M[1])):
        S[j + len(L[1])] = M[1][j]
    #print('Sass',S)
    K = {}
    for i in range(len(L[1]) + len(M[1])):
        for j in range(len(L[1]) + len(M[1])):
            if S[i].intersection(S[j]) != set([]):
                try:
                    K[i].append(j)
                except:
                    K[i] = [j]
    #print('Kayyy',K)

    return J_mult({Af: 1}, dfactor(len(Loops(K,clos_perm)[1]) - 1))

def split_jones(L1,L2,clos_perm,Wr):
    if __name__ == '__main__':
        #print('entered split jones')
        start_time=time.time()
        inputs=[L1,L2]
        #print('ip',inputs)
        with multiprocessing.Pool() as pool:
            #print('entered pool process poly')
            results = pool.map(process_partial_poly, inputs)
            JPOLY = results
            #print('len JPOLY',len(JPOLY))
            #print('L1 L2 inside split done')
            L1_pp=JPOLY[0]
            L2_pp=JPOLY[1]
            #print('L1pp len', len(L1_pp))
            #print('L2pp len', len(L2_pp))
            T=2
            ## GLUING back information from L1_pp and L2_pp. Get bracket polynomial as state sum terms.
            if len(L1_pp)>= len(L2_pp):
                #print('L1>L2')
                chunks = [L1_pp[i:i + len(L1_pp)//T] for i in range(0,len(L1_pp), len(L1_pp)//T)]
                args_list = [(chunk, L2_pp) for chunk in chunks]  
            else:
                #print('L2>L1')
                chunks = [L2_pp[i:i + len(L2_pp)//T] for i in range(0,len(L2_pp), len(L2_pp)//T)]
                args_list = [(chunk, L1_pp) for chunk in chunks]
            #for x in chunks:
            #    print(x)
            #print('-------args_list below-------')
            #for x in args_list:
            #    print(x)
            #print('Argslist len', len(args_list))
            
            tups=pool.map(process_chunk, args_list)
            #print('len tups',len(tups))
            #print(tups)
            #for x in tups:
                #print('tup x',len(x),x[0])
                #print('')
                #for elt in x:
                #    print(elt)
                #print('-----')
            elements_to_process = []
            for x in tups:
                #print('X',x)
                elements_to_process += [elt+(clos_perm,) for elt in x]
            #print('elts to process',elements_to_process)
            #print('elts tp', len(elements_to_process))
            Bkt_poly_summands = pool.map(mult_2_elements, elements_to_process)
            #print('Bkt poly summands', Bkt_poly_summands)
            # Sum up the results from parallel processing
            Bkt_poly={0:0}
            for bpoly in Bkt_poly_summands:
                Bkt_poly=J_add(Bkt_poly,bpoly)
            final_result = J_mult(Bkt_poly,{3*(-Wr):(-1)**abs(Wr)})
            #print('split jones final result:',final_result)
        end_time = time.time()
        time_taken = end_time - start_time
        print('Time taken split parallel :', time_taken)
        return final_result

def jones_execution(num_projj, input_knot):
    '''Uniform Distribution of points on S^2'''
    n = num_projj+2
    points = []
    if n > 1:
        points = fibonacci_sphere(n)[1:-1]
    else:  
        points = [np.array([0.,0.,1.])]
    num_proj=num_projj

    # JPOLY stores Jones Polynomial for each projection.
    JPOLY = []
    # Jnone keeps count of irregular projections
    Jnone = 0
    for proj_vector_i in range(len(points)):
        proj_vector=points[proj_vector_i]
        # Jones polynomial along a specific projection vector
        rj = process_projection((form_indices(input_knot,closed)[0], proj_vector, form_indices(input_knot,closed)[1]))
        if rj!=None:
            JPOLY.append(rj)
        else: 
            # count it as a rejected projection (irregular diagram)
            Jnone += 1
    # Addition of Jones polynomials
    Zpoly = {0:0}
    for i in range(len(JPOLY)):
        Zpoly = J_add(JPOLY[i], Zpoly)
    # Normalize the result by the number of successful Jones polynomial calculations
    print('regular vs. irregular projections', num_proj-Jnone, Jnone)
    num_proj = num_proj - Jnone
    try:
        Zpoly = J_smult(1./num_proj, Zpoly)
    except:
        pass
    # Final Jones polynomial of the input curves
    jones_final= Zpoly
    print('Jones polynomial : ', jones_final)


############### 
''' number of projections'''
num_projj=1
''' input Knot/Link'''
'''ex1 : Trefoil'''
Trefoil=np.array([[[1, 0, 0],[4, 0, 0],[1, 6, 2],[0, 2, -5],[5, 2, 5],[4, 6, -2]]], dtype=object)
'''ex1 : Borromean Rings'''
u1=[[0,0,0],[1,1,0],[2,2,0.5],[3,3,0.5],[4,4,0],[5,5,0],[6,6,0.5],[7,7,0.5],[8,7,0.5],[9,5,0.2],[9,3,0.2],[8,0,0.2],[8,-1,0.2],[6,-1.5,0],[4,-2,0],[2,-1.5,0]]
u2=[[1,0,0.5],[4,0,0],[5,1,0],[5,4,0.5],[4,5,0.5],[3,6,0],[2,7,0],[-1,6,0],[-1,3,0.5]]
u3=[[6,0,0.5],[7,6,0],[6,7,0],[3,7,0.5],[2,6,0.5],[2,3,0],[3,2,0],[4,1,0.5]]
Borr=np.array([u1,u2,u3], dtype=object)
input_knot=Trefoil
## If Closed Curves, then closed = 1. If Open Curves, then closed = 0.
closed=1
## If Parallel Alg, then parallel = 1. If Serial, then parallel= 0.
parallel=1
''' Evaluation of the Jones Polynomial of input_knot'''
jones_execution(num_projj,input_knot)
###############


