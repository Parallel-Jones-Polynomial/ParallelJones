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
import random
from collections import defaultdict
from itertools import combinations
import networkx as nx
from networkx.algorithms.community import kernighan_lin_bisection


''' 
get_jones_poly : 
coords signifies the co-ordinates of points on the curve in 3space
proj_vector signifies a particular direction of projection
inds signifies the connections among the points in the curve in concern.
'''
"""Calculate the Jones polynomial wrt a projection vector."""
def get_jones_poly(coords, proj_vector, Inds, parallel):
    inds=Inds[0]
    clos_perm=Inds[1]
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
    ## Check if the projection is irregular i.e. if 2 or more vertices are overlapping.
    Check_proj = np.unique([str(proj_new[i]) for i in proj_new], return_counts=True)
    if max(Check_proj[1]) >= 2:
        return None
    ## Generate crossing matrices and associated properties.
    bool_mask, over_or_under, right_or_left, inds = crossing_matrices_bool(proj_new, inds, depth)
    ## Calculate Writhe.
    Wr = get_writhe(bool_mask, over_or_under, right_or_left)
    ## Do Reidemeister Moves
    if np.count_nonzero(bool_mask)/2.>=20.:
        bool_mask=simplification(bool_mask, over_or_under,inds)
        Wr = get_writhe(bool_mask, over_or_under, right_or_left)
        if np.count_nonzero(bool_mask)/2.<1:
            return None
        else:
            pass
    ## Splitting Condition
    if parallel==1:
        ## SPLITTING OF THE KNOT INTO 2 LINKOIDS L1 and L2
        def construct_crossing_graph(inds, bool_mask):
            # initialize empty graph
            G = nx.Graph()
            # Collect crossings from bool_mask
            crossings = []
            k=bool_mask.shape[0]
            for i in range(k):
                for j in range(i + 1, k):
                    if bool_mask[i][j]:
                        crossings.append((i, j))
            # Each crossing becomes a vertex in G
            G.add_nodes_from(crossings)
            # add an edge between vertices v1 and v2 if there's a path in the knot diagram b/w v1 and v2
            for c1 in crossings:
                for c2 in crossings:
                    if c1 != c2:
                        if is_connected(c1, c2, inds,crossings):
                            G.add_edge(c1, c2)
            return G

        def is_connected(c1, c2, inds, crossings):
            e1_1 = c1[0]
            e1_2 = c1[1]
            e2_1 = c2[0]
            e2_2 = c2[1]
            # All the edges participating in crossings
            starts=[C[0] for C in crossings]+[C[1] for C in crossings] 
            # Check if e1_start and e2_start are directly connected  
            def traverse_path(start_edge, target_edge):
                current_edge = start_edge     
                conn=False   
                flag=0
                while flag==0:
                    if current_edge in starts:
                        if current_edge == start_edge:
                            current_edge = inds[current_edge]
                        elif current_edge == target_edge:
                            conn= True
                            flag=1
                        else:
                            conn= False
                            flag=1
                    else:
                        # Move to the next edge in the path
                        try:
                            current_edge = inds[current_edge] 
                        except:
                            current_edge = clos_perm[current_edge]    
                return conn
            check_conn = any(traverse_path(x, y) for x in (e1_1, e1_2) for y in (e2_1, e2_2))
            return check_conn

        # Construct the crossing graph
        G = construct_crossing_graph(inds, bool_mask)
        # Apply Kernighan-Lin bisection algorithm
        def get_subgraphs(G, partition):
            subset1, subset2 = partition
            return G.subgraph(subset1).copy(), G.subgraph(subset2).copy()

        def split_ntimes(G, n):
            if n == 0:
                return [G]            
            G1, G2 = get_subgraphs(G,kernighan_lin_bisection(G))
            return split_ntimes(G1, n - 1) + split_ntimes(G2, n - 1)
        
        def create_arguments(n, bool_mask, over_or_under, right_or_left, inds,clos_perm):
            partition = split_ntimes(G, n)
            arcs_per_piece = find_arcs(inds, partition,clos_perm)
            args = []          
            for arcs in arcs_per_piece:
                inds_sub = {k: v for k, v in inds.items() if k in arcs}
                for k in arcs:
                    if k not in inds_sub and k in clos_perm:
                        inds_sub[k] = clos_perm[k]
                modf_inds_sub = {k: [v] for k, v in inds_sub.items()}
                BM_sub = np.copy(bool_mask)
                OU_sub = np.copy(over_or_under)
                RL_sub = np.copy(right_or_left)
                
                for i in range(len(BM_sub)):
                    if i not in arcs:
                        BM_sub[i, :] = False
                        BM_sub[:, i] = False
                        OU_sub[i, :] = False
                        OU_sub[:, i] = False
                        RL_sub[i, :] = False
                        RL_sub[:, i] = False
                
                states = []
                args.append((BM_sub, OU_sub, RL_sub, inds_sub, {}, modf_inds_sub, states, list(np.copy(states))))
            
            return args
        if np.count_nonzero(bool_mask)/2.>7:
            All_args=create_arguments(3, bool_mask, over_or_under, right_or_left, inds, clos_perm) 
        elif 7>=np.count_nonzero(bool_mask)/2.>3:
            All_args=create_arguments(2, bool_mask, over_or_under, right_or_left, inds, clos_perm)
        else: 
            All_args=create_arguments(1, bool_mask, over_or_under, right_or_left, inds, clos_perm)
        ## Return the Jones polynomial in terms of the bracket polynomial states of the L1 and L2 linkoids.
        SJ=split_jones(All_args,clos_perm,Wr)
        return SJ
    else:
        '''Serial Code'''
        ## Compute the Bkt and Jones without splitting.
        states=[]
        modf_inds = {}
        for i in inds:
            modf_inds[i]=[inds[i]]
        start_time=time.time()
        L_bkt = get_partial_poly(bool_mask, over_or_under, right_or_left,inds,clos_perm, modf_inds,states,list(np.copy(states)))
        end_time=time.time()
        return J_mult(L_bkt[0][()],{3*(-Wr):(-1)**abs(Wr)})
        
"""Create crossing matrices, new vertices and other properties """
def crossing_matrices_bool(proj, inds, depth):
    ## Create ordered edges from inds dict
    edges={}
    ct=0
    for i in inds.keys():
        edges[i]=(i,inds[i])
        ct+=1 
    ## Create 4 types of matrices :
    ## bool_mask : a_ij = True (False) if edge i and edge j have (does not have) a crossing. 
    ## over_or_under : a_ij = True (False) if edge i lies above (does not lie above) edge j. 
    ## right_or_left : a_ij = True (False) if cross product: edge i x edge j is positive (negative). (Right hand thumb rule).
    ## u :  a_ij = coords of intersection point (if any) b/w edge i and edge j.
    modf_inds0 = {}
    for i in inds:
        modf_inds0[i]=[inds[i]]
    Loops_as_sets=Loops(modf_inds0,{},0)[1]
    total_inds=[]
    for elt in Loops_as_sets:
        total_inds+=list(elt)
    I_inds=len(total_inds)
    bool_mask = np.zeros((I_inds,I_inds), dtype=bool)
    over_or_under = np.zeros((I_inds,I_inds), dtype=bool)
    right_or_left = np.zeros((I_inds,I_inds), dtype=bool)
    u = np.zeros((I_inds,I_inds), dtype=object)
    ## dif : stores vector rep for any edge, (p1,p2), in terms of coords of the endpoints : vec(p2)-vec(p1).
    ## flip neg : similar to dif but replace [x,y] -> [y,-x] (flipped and negative). Used in Cramer rule later.
    dif={}
    flip_neg={}
    for e in edges.keys():
        dif[e]=[proj[edges[e][1]][0]-proj[edges[e][0]][0],proj[edges[e][1]][1]-proj[edges[e][0]][1]]
        flip_neg[e]=[proj[edges[e][1]][1]-proj[edges[e][0]][1],-proj[edges[e][1]][0]+proj[edges[e][0]][0]]
    ## Iterate through each pair of edges to find intersection point(s) and update the matrices.
    ## intersection pt vs edges which contain it.
    intp_vs_edges = {} 
    num_edges=edges.keys()
    for ii in range(len(num_edges)):
        for jj in range( ii+1,len(num_edges)):
            i=list(num_edges)[ii]
            j=list(num_edges)[jj]
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
            ## If an intersection is found, matrices are updated as follows: 
            if int_pt != None:
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
                ## Check if intersection point is a crossing i.e. lies within the edge segments.
                if zero_val<t<1-zero_val and zero_val<s<1-zero_val: 
                    bool_mask[i,j]=True
                    bool_mask[j,i]=True
                    intp_vs_edges[str(int_pt)]=[i,j]
                    ## Record which edge is over the other at the intersection point in over_or_under matrix.
                    zi=depth[edges[i][0]]+t*(depth[edges[i][1]]-depth[edges[i][0]])
                    zj=depth[edges[j][0]]+s*(depth[edges[j][1]]-depth[edges[j][0]])
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
    ## Update the matrices to account for the newly introduced vertices.
    ## BM_2, OU_2 and RL_2 are expanded to make space for subdivided edges.
    modf_inds1 = {}
    for i in inds:
        modf_inds1[i]=[inds[i]]
    Loops_as_sets1=Loops(modf_inds1,{},0)[1]
    total_inds1=[]
    for elt in Loops_as_sets1:
        total_inds1+=list(elt)
    J_inds=len(total_inds1)
    BM_2= np.zeros((J_inds,J_inds), dtype=bool)
    OU_2= np.zeros((J_inds,J_inds), dtype=bool)
    RL_2= np.zeros((J_inds,J_inds), dtype=bool)
    for pt in intp_vs_newedges:
        i=intp_vs_newedges[pt][0]
        j=intp_vs_newedges[pt][1]
        BM_2[i,j]=True
        BM_2[j,i]=True
        i0=intp_vs_edges[pt][0]
        j0=intp_vs_edges[pt][1]
        OU_2[i,j]=np.copy(over_or_under)[i0,j0]
        OU_2[j,i]=np.copy(over_or_under)[j0,i0]
        RL_2[i,j]=np.copy(right_or_left)[i0,j0]
        RL_2[j,i]=np.copy(right_or_left)[j0,i0] 
    return BM_2, OU_2, RL_2, inds

def gauss_code(bool_mask, over_or_under,inds):
    Gc=[]
    labelss={}
    sl=0
    k=bool_mask.shape[0]
    binss=[]
    i=0
    while i not in binss:
        for j in range(k):
            if bool_mask[i,j]==True:
                if '{}{}'.format(i,j) not in labelss:
                    sl=sl+1
                    labelss['{}{}'.format(i,j)]=sl
                    labelss['{}{}'.format(j,i)]=sl
                if over_or_under[i,j]==True:
                    Gc.append(labelss['{}{}'.format(i,j)])
                else:
                    Gc.append(-labelss['{}{}'.format(i,j)])
        binss.append(i)
        i=inds[i]
    return Gc

"""Calculate BRACKET POLYNOMIAL using crossing matrices and associated properties"""
def get_partial_poly(bool_mask, over_or_under, right_or_left, inds, clos_perm, modf_inds,states,S):
    ## Check any entry in bool_mask is True
    if np.any(bool_mask):
        edge1 = np.argmax(np.any(bool_mask, 0))
        edge2 = np.argmax(bool_mask[edge1,:])
        
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
        if over_or_under[edge1, edge2] == right_or_left[edge1, edge2]:
            modf_inds1[edge1].append(inds[edge2])
            modf_inds1[edge2].append(inds[edge1])
            modf_inds2[edge1].append(edge2)
            try:
                modf_inds2[inds[edge1]].append(inds[edge2]) 
            except:
                modf_inds2[inds[edge1]]=[inds[edge2]]
        else:
            modf_inds1[edge1].append(edge2) 
            try:
                modf_inds1[inds[edge1]].append(inds[edge2])
            except:
                modf_inds1[inds[edge1]]=[inds[edge2]]
            
            modf_inds2[edge1].append(inds[edge2])
            modf_inds2[edge2].append(inds[edge1])
        states1.append(1)
        states2.append(-1)
        get_partial_poly(bool_mask1, np.copy(over_or_under), np.copy(right_or_left), inds, clos_perm, modf_inds1,states1,S)
        get_partial_poly(bool_mask2, np.copy(over_or_under), np.copy(right_or_left), inds, clos_perm, modf_inds2,states2,S)
    else:
        raw_summand=(sum(states),Loops(modf_inds,clos_perm,1)[1])
        if not S:
            S.append({})
        key = tuple(sorted(map(lambda s: tuple(sorted(s)), (x for x in raw_summand[1] if x)))) #tuple(sorted(map(tuple, (s for s in raw_summand[1] if s))))
        # Count the number of empty sets
        empty_count = sum(1 for s in raw_summand[1] if not s)
        # Create the entry value
        entry_value = (raw_summand[0], empty_count)       
        # Update the dictionary in A[0]
        if key==():
            if key not in S[0]:
                S[0][key] = J_mult({entry_value[0]:1},dfactor(entry_value[1]-1)) 
            else:
                S[0][key] = J_add(S[0][key],J_mult({entry_value[0]:1},dfactor(entry_value[1]-1)))
        else:
            if key not in S[0]:
                S[0][key] = J_mult({entry_value[0]:1},dfactor(entry_value[1])) 
            else:
                S[0][key] = J_add(S[0][key],J_mult({entry_value[0]:1},dfactor(entry_value[1])))
    return S

''' WRITHE '''
def get_writhe(bool_mask, over_or_under, right_or_left):
    BM=np.copy(bool_mask)
    Wr=(np.sum(np.int32(over_or_under == right_or_left)[BM])-np.sum(np.int32(over_or_under != right_or_left)[BM]))/2.
    return Wr

# Find the arc that forms the backbone of L1.
## I is the backbone arc and J stores the auxillary arcs.
def find_arcs(inds, partition,clos_perm):
    subgraph_edges = []
    for subgraph in partition:
        C_bag = list(subgraph.nodes())
        I = [C[0] for C in C_bag] + [C[1] for C in C_bag]
        I_final = I[:]
        I_complement = [node[0] for other_subgraph in partition if other_subgraph != subgraph for node in other_subgraph.nodes()]+[node[1] for other_subgraph in partition if other_subgraph != subgraph for node in other_subgraph.nodes()]
        for i in I:
            print('start i', i)
            flag = 0
            visited = set()  # Track visited values to detect cycles
            
            while flag == 0:
                if i in visited:  # Cycle detected, break to prevent infinite loop
                    print(f"Cycle detected at {i}, breaking loop.")
                    break
                visited.add(i)  # Mark as visited
                
                try:
                    i = inds[i]  # Try updating i using inds
                    print('try 1', i)
                except KeyError:
                    try:
                        i = clos_perm[i]  # Try updating i using clos_perm
                        print('try 2', i)
                    except KeyError:
                        i += 1  # If both fail, increment i
                        print('try 3', i)

                if i not in I_complement:
                    if i not in I_final:
                        I_final.append(i)
                else:
                    flag = 1  # Stop if i is in I_complement
        print('I_final after', I_final)
        subgraph_edges.append(I_final)
        print('final edges', subgraph.edges())
    return subgraph_edges


## Designates unique index numbers to each coordinate of the polygonal knot.
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
    # knot_data : indices dictionary, closure permutations, particular projection vector
    #print('KD',knot_data)
    rj = get_jones_poly(*knot_data)
    #print('rj',rj)
    return rj

def process_partial_poly(knot_data):
    #print('entered process partial poly')
    rj = get_partial_poly(*knot_data)
    #print('rj process partial poly:', rj)
    return rj


# ALL THAT IS TO BE PICKLED NEED TO BE DEFINED OUTSIDE
def process_chunk(args):
    """Generates all possible pairs from W1_chunk and W2"""
    W1_chunk, W2 = args
    W1_items = W1_chunk.items()  
    W2_items = W2.items()
    R = [(w1_item, w2_item) for w1_item in W1_items for w2_item in W2_items]
    return R

def f(t):
    '''Return transpositions keeping only unpaired vertices by absorbing paired vertices'''
    flat_list = [item for sublist in t for item in sublist]
    # Count the occurrences of each element
    count = defaultdict(int)
    for num in flat_list:
        count[num] += 1
    # Filter the list to include only numbers that appear once
    return tuple(num for num in flat_list if count[num] == 1)

def mult_2_elements(L_M):
    '''Multiplication of 2 states'''
    L, M, clos_perm, flag = L_M
    vertices = L[0] + M[0]
    wgt_factor=J_mult(L[1],M[1])
    G = nx.Graph()
    # Add nodes with a unique identifier (Since same node may appear in both pieces and we want to treat them separately)
    node_map = [(v, i) for i, v in enumerate(vertices)]
    G.add_nodes_from(node_map)  # Add multiple nodes at once
    # Add edges between vertices that share common points
    node_sets = {idx: set(v) for v, idx in node_map}
    # Generate all edges efficiently
    edges = [((v1, idx1), (v2, idx2)) for (v1, idx1), (v2, idx2) in combinations(node_map, 2) if node_sets[idx1] & node_sets[idx2]]
    G.add_edges_from(edges)
    # Find connected components using networkx
    components = list(nx.connected_components(G))
    # Apply the function f(t) to each component and construct the collection of transposition P
    P = tuple(f([node[0] for node in component]) for component in components)
    #print('Resultant', P)
    null_c = sum(1 for t in P if not t)
    P1 = tuple(t for t in P if t)
    if flag==0:
        resultant={P1 : J_mult(wgt_factor, dfactor(null_c))}
    else:
        # When Mult is done in the last step. d^{-1} factor has to be accounted for
        resultant={P1 : J_mult(wgt_factor, dfactor(null_c-1))}
    #print('FINAL RESULTANT', resultant)
    return resultant

def chunk_list(data, T):
    """Splits a dictionary into T approximately equal chunks."""
    items = list(data.items())
    chunk_size = max(1, len(items) // T)  # Ensure chunk size is at least 1
    return [dict(items[i:i + chunk_size]) for i in range(0, len(items), chunk_size)]


def split_jones(All_args,clos_perm,Wr):
    '''Splitting, Gluing and Parallelisation CODE BLOCK'''
    if __name__ == '__main__':
        inputs=All_args
        with multiprocessing.Pool() as pool:
            '''Parallel-1 : OBTAIN STATE SPACES FOR EACH LINKOID PIECE IN PARALLEL'''
            results = pool.map(process_partial_poly, inputs)
            results_array = np.array(results)
            L_pp = results_array[:len(All_args), 0]
            print('LPP',L_pp)
            ''' GLUING OF STATES IN PARALLEL''' 
            def chunk_list(data, T):
                """Splits a dictionary into T approximately equal chunks."""
                items = list(data.items())
                chunk_size = max(1, len(items) // T)  # Ensure chunk size is at least 1
                return [dict(items[i:i + chunk_size]) for i in range(0, len(items), chunk_size)]

            def multiply_lists(list1, list2, T, flag):
                """Performs parallel chunk-wise multiplication of list1 and list2."""               
                # Ensure list1 is always the larger one by swapping to maintain consistency
                if len(list1) < len(list2):
                    list1, list2 = list2, list1
                # Split larger list into T chunks
                chunks = chunk_list(list1, T)               
                # Prepare arguments for parallel processing
                args_list = [(chunk, list2) for chunk in chunks]
                '''Parallel-2 : Parallel creation of pairs from all the chunks'''
                tups = pool.map(process_chunk, args_list)
                elements_to_process = [elt + (clos_perm,) + (flag,) for x in tups for elt in x]
                '''Parallel-3 : Parallel creation of pairs from all the chunks'''
                return pool.map(mult_2_elements, elements_to_process)
            # MULTIPLICATION TAKES PLACE IN LEVELS AS FOLLOWS
            Num_levels=math.log(len(All_args))/math.log(2)
            print('Num Levels', Num_levels)
            level_list=L_pp
            #print('INITIAL Level List', level_list)
            T=2
            for i in range(int(Num_levels)): 
                print('Level Number',i)
                ultimate_level_list=[]
                # flag1 = 1 only if we are at the FINAL LEVEL
                flag1 = 0 if i != int(Num_levels) - 1 else 1
                # Mult of two blocks in a given level
                level_L = [multiply_lists(level_list[2*j], level_list[2*j+1], T, flag1) for j in range(len(level_list) // 2)]
                level_list=level_L
                print('Next Level List', level_list)
                for pplist in level_list:
                    temp_dict = {}
                    for dict_elt in map(dict, pplist): 
                        key1, value = next(iter(dict_elt.items()))
                        key = tuple(sorted(tuple(sorted(s)) for s in key1 if s))
                        temp_dict[key] = J_add(temp_dict[key], value) if key in temp_dict else value
                    ultimate_level_list.append(temp_dict)
                level_list=ultimate_level_list
                T=2*T
                print('Level List at this level', ultimate_level_list)
            try:    
                final_result=J_mult(level_list[0][()],{3*(-Wr):(-1)**abs(Wr)})
            except:
                def count_loops(clos_perm, target_tuple):
                    G = nx.Graph()
                    G.add_edges_from(target_tuple) 
                    # Apply clos_perm to update connections
                    for a, b in clos_perm.items():
                        if G.has_node(a):
                            G = nx.relabel_nodes(G, {a: b})
                    # Count the number of connected components (loops)
                    return nx.number_connected_components(G)
                final_result={0:0}
                for tupp in level_list[0]:
                    final_result=J_add(final_result,J_mult(level_list[0][tupp],dfactor(int(count_loops(clos_perm,tupp)))))
                final_result=J_mult(final_result,{3*(-Wr):(-1)**abs(Wr)})
            print('split jones final result:',final_result)
        return final_result
    
def jones_execution(num_projj, input_knot, closed, parallel):
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
        rj = process_projection((form_indices(input_knot,closed)[0], proj_vector, form_indices(input_knot,closed)[1],parallel))
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



import argparse
import ast

def str_to_np_array(s):
    try:
        return np.array(ast.literal_eval(s), dtype=object)
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Invalid array format: {e}")

parser = argparse.ArgumentParser(description="Run jones_execution with command-line arguments.")

parser.add_argument("num_projj", type=int, help="Number of projections")
parser.add_argument("input_knot", type=str, help="Input knot as a string (e.g., '[[[1, 0, 0],[4, 0, 0],...]]')")
parser.add_argument("closed", type=int, choices=[0, 1], help="1 if knot is closed, 0 otherwise")
parser.add_argument("parallel", type=int, choices=[0, 1], help="1 for parallel execution, 0 otherwise")

args = parser.parse_args()

# Convert input_knot from string to NumPy array
input_knot_array = str_to_np_array(args.input_knot)

# Convert closed and parallel to bool
closed_bool = bool(args.closed)
parallel_bool = bool(args.parallel)

# Run the function
result = jones_execution(args.num_projj, input_knot_array, closed_bool, parallel_bool)

# Output the result
print(result)
