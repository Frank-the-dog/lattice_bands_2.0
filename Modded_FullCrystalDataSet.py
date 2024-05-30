#Full Crystal Data Set

#Input_Crystal of the style "FCC" or "cF" (is case sensative!)
# example: dim = (2,3,5,0.2,0.7,0.8)
# The dimentions (a,b,c,alpha,beta,gamma). Only those defined in a particular crystal structure will be used be the program. Angles are in Radians.

# Crystal definitions taken from https://arxiv.org/pdf/1004.2974.pdf

import numpy as np

CrystalList = np.array(\
[                               #The source data for the lattice selection of form [text,name,set dim dimentions,mod number]
                                #Mod number is for called a couple lattices that cannot be gotted with an inequality
    ["Cubic (CUB,cP)","CUB",0.8,0.5,01.0,0.6,0.4,0.9," "],
    ["Face-centered Cubic (FCC,cF)","FCC",0.8,0.5,01.0,0.6,0.4,0.9," "],
    ["Body-centered Cubic (BCC,cI)","BCC",0.8,0.5,01.0,0.6,0.4,0.9," "],
    ["Tetragonal (TET,tP)","TET",0.8,0.5,01.2,0.6,0.4,0.9," "],#
    ["Body-centered Tetragonal form 1 (BCT_1,tI)","BCT",1.0,0.5,0.8,0.6,0.4,0.9," "],
    ["Body-centered Tetragonal form 2 (BCT_2,tI)","BCT",0.8,0.5,01.0,0.6,0.4,0.9," "],
    ["Orthorhombic (ORC,oP)","ORC",0.8,0.5,01.0,0.6,0.4,0.9," "],#
    ["Face-centered Orthorhombic form 1 (ORCF_1,oF)","ORCF",0.5,0.8,01.0,0.6,0.4,0.9," "],#
    ["Face-centered Orthorhombic form 2 (ORCF_2,oF)","ORCF",0.7,0.8,0.9,0.6,0.4,0.9," "],#
    ["Face-centered Orthorhombic form 3 (ORCF_3,oF)","ORCF",0.5979,0.8,0.9,0.6,0.4,0.9,"3"],#
    ["Body-centered Orthorhombic (ORCI,oI)","ORCI",0.8,0.5,01.0,0.6,0.4,0.9," "],#
    ["C-centered Orthorhombic (ORCC,oS)","ORCC",0.5,0.8,01.0,0.6,0.4,0.9," "],#
    ["Hexagonal (HEX,hP)","HEX",0.8,0.5,01.0,0.6,0.4,0.9," "],#
    ["Rhombohedral form 1 (RHL_1,hR)","RHL",0.5,0.6,01.0,1.9,0.4,0.9," "],#
    ["Rhombohedral form 2 (RHL_2,hR)","RHL",0.5,0.6,01.0,1.2,0.4,0.9," "],#
    ["Monoclinic (MCL,mP)","MCL",0.9,0.8,01.0,1.0,0.4,0.9," "],#
    ["C-centered Monoclinic form 1 (MCLC_1,mS)","MCLC",0.4,0.7,01.0,1.2,0.4,0.9," "],#
    ["C-centered Monoclinic form 2 (MCLC_2,mS)","MCLC",0.5,0.6,01.0,0.9,0.4,0.9,"2"],#
    ["C-centered Monoclinic form 3 (MCLC_3,mS)","MCLC",0.5,0.6,01.0,0.9,0.4,0.9," "],#
    ["C-centered Monoclinic form 4 (MCLC_4,mS)","MCLC",0.9,0.8,01.0,1.0,0.4,0.9,"4"],#
    ["C-centered Monoclinic form 5 (MCLC_5,mS)","MCLC",0.9,0.8,01.0,0.9,0.4,0.9," "],#
    ["Triclinic forms 1a,1b,2a,2b (TRI,aP)","TRI",0.9,0.8,01.0,1.0,0.4,0.9," "],#
])



def Fullcrystaldata(dim,Input_Crystal,modnumber):
    b1 = (None,None,None)

    # So this part looks at the Input_Crystal id, and have a calculation section for each crystal
    # All the data is literally pulled straight from https://arxiv.org/pdf/1004.2974.pdf

    # I'll comment on the first one so you know what everything does
    if Input_Crystal == "CUB" or Input_Crystal == "cP": 
        a1 = (dim[0],0,0)    #Finds primative lattice vectors from standard lattice vectors
        a2 = (0,dim[0],0)
        a3 = (0,0,dim[0])
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))   #Finds momentum space lattice vectors from primative lattice vectors
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        Corners = np.zeros((4,3))
        Corners[0] = 0.0*b1+0.0*b2+0.0*b3    # Calculate the point positions based on the momentum lattice vectors
        Corners[1] = 0.5*b1+0.5*b2+0.0*b3    # The factors are taken straight from the source paper
        Corners[2] = 0.5*b1+0.5*b2+0.5*b3
        Corners[3] = 0.0*b1+0.5*b2+0.0*b3
        Segments = np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])  #Which points to connect with plot lines
        Plot_List = [0,3,1,0,2,3,"break",1,2]           #What path the band plot follows
        Corner_Labels = ["$\Gamma$","M","R","X"]        #The letter names for each corner. $\Gamma$ calls laTex code, and works fine, but generates an error.
        Ghosted_Segments = []                   #The dashed (non-sampled) line segments
        Outer_Segments = [[1,2]]                #The line segments that are part of the outer (full Brillouin zone) frame
        unfold_instructions = [["mirror",[0,1,2]],["mirror",[0,2,3]],["mirror",[0,1,3]],["2fmirror",[[0,1,2],[0,2,3]]],\
                               ["2fmirror",[[0,1,3],[0,1,2]]],["2fmirror",[[0,2,3],[0,1,3]]],["2fmirror",[[0,1,3],[0,1,2]]],
                               ["mirror",[0,1,3]],["mirror",[0,2,3]]]       #The mirror commands to get the outer segments to build a full frame

    if Input_Crystal == "FCC" or Input_Crystal == "cF":
        a1 = (0,0.5*dim[0],0.5*dim[0])
        a2 = (0.5*dim[0],0,0.5*dim[0])
        a3 = (0.5*dim[0],0.5*dim[0],0)
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        Corners = np.zeros((6,3))
        Corners[0] = 0.0*b1+0.0*b2+0.0*b3
        Corners[1] = 3/8*b1+3/8*b2+3/4*b3
        Corners[2] = 0.5*b1+0.5*b2+0.5*b3
        Corners[3] = 5/8*b1+1/4*b2+5/8*b3
        Corners[4] = 1/2*b1+1/4*b2+3/4*b3
        Corners[5] = 1/2*b1+0.0*b2+1/2*b3
        Segments = np.array([[0,1],[0,2],[0,5],[1,2],[1,4],[2,3],[2,4],[3,4],[3,5],[4,5]])
        Plot_List = [0,5,4,1,0,2,3,4,2,1,'break',3,5]
        Corner_Labels = ['$\Gamma$','K','L','U','W','X']
        Ghosted_Segments = []
        Outer_Segments = [[1,4],[3,4]]
        unfold_instructions = [["mirror",[0,1,4]],["mirror",[0,2,3]],["mirror",[0,1,2]],["2fmirror",[[0,1,4],[0,2,3]]],\
                               ["2fmirror",[[0,1,4],[0,1,2]]],["2fmirror",[[0,2,3],[0,1,2]]],["2fmirror",[[0,2,3],[0,1,4]]],]

    if Input_Crystal == "BCC" or Input_Crystal == "cI":
        a1 = (-0.5*dim[0],0.5*dim[0],0.5*dim[0])
        a2 = (0.5*dim[0],-0.5*dim[0],0.5*dim[0])
        a3 = (0.5*dim[0],0.5*dim[0],-0.5*dim[0])
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        Corners = np.zeros((4,3))
        Corners[0] = 0.0*b1+0.0*b2+0.0*b3
        Corners[1] = 1/2*b1-1/2*b2+1/2*b3
        Corners[2] = 1/4*b1+1/4*b2+1/4*b3
        Corners[3] = 0.0*b1+0.0*b2+1/2*b3
        Segments = np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
        Plot_List = [0,1,3,0,2,1,'break',2,3]
        Corner_Labels = ['$\Gamma$','H','P','N']
        Ghosted_Segments = []
        Outer_Segments = [[1,2]]
        unfold_instructions = [["mirror",[0,2,3]],["mirror",[0,3,1]],["mirror",[0,2,1]],["2fmirror",[[0,1,2],[0,2,3]]],\
                               ["2fmirror",[[0,2,3],[0,1,2]]],["2fmirror",[[0,1,3],[0,1,2]]],["2fmirror",[[0,2,3],[0,1,3]]],
                               ["mirror",[0,2,3]]]

    if Input_Crystal=="TET" or Input_Crystal=="tP":
        a1 = (dim[0],0,0)
        a2 = (0,dim[0],0)
        a3 = (0,0,dim[2])
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        Corners = np.zeros((6,3))
        Corners[0] = 0.0*b1+0.0*b2+0.0*b3
        Corners[1] = 1/2*b1+1/2*b2+1/2*b3
        Corners[2] = 1/2*b1+1/2*b2+0.0*b3
        Corners[3] = 0.0*b1+1/2*b2+1/2*b3
        Corners[4] = 0.0*b1+1/2*b2+0.0*b3
        Corners[5] = 0.0*b1+0.0*b2+1/2*b3
        Segments = np.array([[0,2],[0,4],[0,5],[1,2],[1,3],[1,5],[2,4],[3,4],[3,5]])
        Plot_List = [0,4,2,0,5,3,1,5,'break',4,3,'break',2,1]
        Corner_Labels = ['$\Gamma$','A','M','R','X','Z']
        Ghosted_Segments = []
        Outer_Segments = [[1,2],[1,3]]
        unfold_instructions = [["mirror",[0,5,1]],["mirror",[0,5,3]],["mirror",[0,2,4]],["2fmirror",[[0,1,5],[0,5,3]]]]

    if Input_Crystal=="BCT" or Input_Crystal=="tI":
        a1 = (-0.5*dim[0],0.5*dim[0],0.5*dim[2])
        a2 = (0.5*dim[0],-0.5*dim[0],0.5*dim[2])
        a3 = (0.5*dim[0],0.5*dim[0],-0.5*dim[2])
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        if dim[2] < dim[0]:    #BCT1
            eta = (1 + dim[2]**2/dim[0]**2)/4
            Corners = np.zeros((7,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = -1/2*b1+1/2*b2+1/2*b3
            Corners[2] = 0.0*b1+1/2*b2+0.0*b3
            Corners[3] = 1/4*b1+1/4*b2+1/4*b3
            Corners[4] = 0.0*b1+0.0*b2+1/2*b3
            Corners[5] = eta*b1+eta*b2-eta*b3
            Corners[6] = -eta*b1+(1-eta)*b2+eta*b3
            Segments = np.array([[0,1],[0,4],[0,5],[1,4],[1,6],[2,3],[2,6],[3,4],[3,5]])
            Plot_List = [0,4,1,0,5,3,2,6,1,"break",4,3]
            Corner_Labels = ["$\Gamma$",'M','N','P','X','Z','$Z_1$']
            Ghosted_Segments = [[2,5],[3,6]]
            Outer_Segments = [[1,6],[3,6],[3,5]]
            unfold_instructions = [["mirror",[1,0,4]],["mirror",[0,3,4]],["mirror",[1,0,6]],["2fmirror",[[0,3,4],[0,1,2]]]]
        else:       #BCT2
            eta = (1 + dim[0]**2/dim[2]**2)/4
            zeta = dim[0]**2/2/(dim[2]**2)
            Corners = np.zeros((9,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = 0.0*b1+1/2*b2+0.0*b3
            Corners[2] = 1/4*b1+1/4*b2+1/4*b3
            Corners[3] = -eta*b1+eta*b2+eta*b3
            Corners[4] = eta*b1+(1-eta)*b2-eta*b3
            Corners[5] = 0.0*b1+0.0*b2+1/2*b3
            Corners[6] = -zeta*b1+zeta*b2+1/2*b3
            Corners[7] = 1/2*b1+1/2*b2-zeta*b3
            Corners[8] = 1/2*b1+1/2*b2-1/2*b3
            Segments = np.array([[0,3],[0,5],[0,8],[1,2],[1,4],[2,5],[2,7],[3,6],[4,8],[5,6],[7,8]])
            Plot_List = [0,5,6,3,0,8,4,1,2,7,8,"break",5,2]
            Corner_Labels = ["$\Gamma$",'N','P','$\Sigma$','$\Sigma_1$','X','Y','$Y_1$','Z']
            Ghosted_Segments = [[1,3],[2,6],[4,7]]
            Outer_Segments = [[2,6],[2,7],[3,6],[4,7]]
            unfold_instructions = [["mirror",[0,2,5]],["mirror",[0,3,6]],["mirror",[0,3,1]],["2fmirror",[[0,5,2],[0,3,4]]]]


    if Input_Crystal == "ORC" or Input_Crystal == "oP":
        a1 = (dim[0],0,0)
        a2 = (0,dim[1],0)
        a3 = (0,0,dim[2])
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        Corners = np.zeros((8,3))
        Corners[0] = 0.0*b1+0.0*b2+0.0*b3
        Corners[1] = 0.5*b1+0.5*b2+0.5*b3
        Corners[2] = 0.5*b1+0.5*b2+0.0*b3
        Corners[3] = 0.0*b1+0.5*b2+0.5*b3
        Corners[4] = 1/2*b1+0.0*b2+1/2*b3
        Corners[5] = 1/2*b1+0.0*b2+0.0*b3
        Corners[6] = 0.0*b1+0.5*b2+0.0*b3
        Corners[7] = 0.0*b1+0.0*b2+1/2*b3
        Segments = np.array([[0,5],[0,6],[0,7],[1,2],[1,3],[1,4],[2,5],[2,6],[3,6],[3,7],[4,5],[4,7]])
        Plot_List =[0,5,2,6,0,7,4,1,3,7,"break",6,3,"break",4,5,"break",2,1]
        Corner_Labels = ['$\Gamma$','R','S','T','U','X','Y','Z']
        Ghosted_Segments = []
        Outer_Segments = [[4,1],[1,2],[1,3]]
        unfold_instructions = [["mirror",[0,4,5]],["mirror",[0,2,6]],["mirror",[0,3,6]]]
    

    if Input_Crystal=="ORCF" or Input_Crystal=="oF":      
        [a,b,c,alpha,beta,gamma] = dim
        a1 = (0, b/2, c/2)
        a2 = (a/2, 0, c/2)
        a3 = (a/2, b/2, 0)
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))

        if modnumber == "3":    #ORCF3    (Never actually triggers)
            zeta = (1+a**2/b**2-a**2/c**2)/4
            eta = (1 + a**2/b**2+a**2/c**2)/4
            Corners = np.zeros((8,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = 1/2*b1+(1/2+zeta)*b2+zeta*b3
            Corners[2] = 1/2*b1+(1/2-zeta)*b2+(1-zeta)*b3
            Corners[3] = 1/2*b1+1/2*b2+1/2*b3
            Corners[4] = 1.0*b1+1/2*b2+1/2*b3
            Corners[5] = 0.0*b1+eta*b2+eta*b3
            Corners[6] = 1/2*b1+0*b2+1/2*b3
            Corners[7] = 1/2*b1+1/2*b2+0*b3
            Segments = np.array([[0,3],[0,5],[0,6],[0,7],[1,5],[1,7],[2,5],[2,6],[4,6],[4,7]])
            Plot_List = [0,6,4,7,0,5,2,6,"break",5,1,7,"break",3,0]
            Corner_Labels = ["$\Gamma$",'A','$A_1$','L','T','X','Y','Z']
            Ghosted_Segments = [[1,4],[2,4]]
            Outer_Segments = [[1,4],[2,4],[5,2],[5,1]]
            unfold_instructions = [["mirror",[1,0,5]],["mirror",[0,7,4]],["mirror",[5,0,2]]]

        elif 1/(dim[0]**2) > (1/(dim[1]**2) + 1/(dim[2]**2)):    #ORCF1 
            zeta = (1+a**2/b**2-a**2/c**2)/4
            eta = (1 + a**2/b**2+a**2/c**2)/4
            Corners = np.zeros((9,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = 1/2*b1+(1/2+zeta)*b2+zeta*b3
            Corners[2] = 1/2*b1+(1/2-zeta)*b2+(1-zeta)*b3
            Corners[3] = 1/2*b1+1/2*b2+1/2*b3
            Corners[4] = 1.0*b1+1/2*b2+1/2*b3
            Corners[5] = 0.0*b1+eta*b2+eta*b3
            Corners[6] = 1*b1+(1-eta)*b2+(1-eta)*b3
            Corners[7] = 1/2*b1+0*b2+1/2*b3
            Corners[8] = 1/2*b1+1/2*b2+0*b3
            Segments = np.array([[0,3],[0,5],[0,7],[0,8],[1,5],[1,8],[2,5],[2,7],[4,6],[4,7],[4,8]])
            Plot_List = [0,7,4,8,0,5,2,7,"break",4,6,"break",5,1,8,"break",3,0]
            Corner_Labels = ["$\Gamma$",'A','$A_1$','L','T','X','$X_1$','Y','Z']
            Ghosted_Segments = [[1,6],[2,6]]
            Outer_Segments = [[1,6],[2,6],[4,6],[1,5],[2,5]]
            unfold_instructions = [["mirror",[0,8,5]],["mirror",[0,2,7]],["mirror",[0,4,7]]]

        else:       #ORCF2
            eta = (1 + a**2/b**2 - a**2/c**2)/4
            phi = (1 + c**2/b**2 - c**2/a**2)/4
            delta = (1 + b**2/a**2 - b**2/c**2)/4
            Corners = np.zeros((11,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = 1/2*b1+(1/2-eta)*b2+(1-eta)*b3
            Corners[2] = 1/2*b1+(1/2+eta)*b2+eta*b3
            Corners[3] = (1/2-delta)*b1+1/2*b2+(1-delta)*b3
            Corners[4] = (1/2+delta)*b1+1/2*b2+delta*b3
            Corners[5] = 1/2*b1+1/2*b2+1/2*b3
            Corners[6] = (1-phi)*b1+(1/2-phi)*b2+1/2*b3
            Corners[7] = phi*b1+(1/2+phi)*b2+1/2*b3
            Corners[8] = 0*b1+1/2*b2+1/2*b3
            Corners[9] = 1/2*b1+0*b2+1/2*b3
            Corners[10] = 1/2*b1+1/2*b2+0*b3
            Segments = np.array([[0,5],[0,8],[0,9],[0,10],[1,3],[1,6],[1,9],[2,10],[3,8],[4,6],[4,10],[6,9],[7,8]])
            Plot_List = [0,9,1,3,8,0,10,4,6,1,"break",2,10,"break",8,7,"break",6,9,"break",5,0]
            Corner_Labels = ["$\Gamma$",'C','$C_1$','D','$D_1$','L','H','$H_1$','X','Y','Z']
            Ghosted_Segments = [[2,7],[3,7],[2,4]]
            Outer_Segments = [[1,6],[6,4],[4,2],[2,7],[7,3],[3,1]]
            unfold_instructions = [["mirror",[0,7,2]],["mirror",[0,1,3]],["mirror",[0,4,6]]]



    if Input_Crystal == "ORCI" or Input_Crystal == "oI":   
        [a,b,c,alpha,beta,gamma] = dim
        a1 = (-a/2, b/2, c/2)
        a2 = (a/2, -b/2, c/2)
        a3 = (a/2, b/2, -c/2)
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        zeta = (1+a**2/c**2)/4
        eta = (1+b**2/c**2)/4
        delta = (b**2-a**2)/(4*c**2)
        mu = (a**2+b**2)/(4*c**2)
        Corners = np.zeros((13,3))
        Corners[0] = 0.0*b1+0.0*b2+0.0*b3
        Corners[1] = -mu*b1+mu*b2+(1/2-delta)*b3
        Corners[2] = mu*b1-mu*b2+(1/2+delta)*b3
        Corners[3] = (1/2-delta)*b1+(1/2+delta)*b2-mu*b3
        Corners[4] = 0.0*b1+1/2*b2+0.0*b3
        Corners[5] = 1/2*b1+0.0*b2+0*b3
        Corners[6] = 0.0*b1+0.0*b2+1/2*b3
        Corners[7] = 1/4*b1+1/4*b2+1/4*b3
        Corners[8] = -zeta*b1+zeta*b2+zeta*b3
        Corners[9] = zeta*b1+(1-zeta)*b2-zeta*b3
        Corners[10] = eta*b1-eta*b2+eta*b3
        Corners[11] = (1-eta)*b1+eta*b2-eta*b3
        Corners[12] = 1/2*b1+1/2*b2-1/2*b3
        Segments = np.array([[0,8],[0,10],[0,12],[1,6],[1,8],[2,10],[4,7],[4,9],[5,7],[5,10],[6,7],[9,12],[11,12]])
        Plot_List = [0,8,1,6,7,4,9,12,0,10,5,7,"break",2,10,"break",11,12]
        Corner_Labels = ['$\Gamma$','L','$L_1$','$L_2$','R','S','T','W','X','$X_1$','Y','Y_1','Z']
        Ghosted_Segments = [[4,8],[1,7],[2,6],[2,7],[5,11],[11,3],[3,9],[3,7]]
        Outer_Segments = [[9,3],[3,11],[3,7],[7,1],[7,2],[2,10],[8,1]]
        unfold_instructions = [["mirror",[0,8,4]],["mirror",[0,1,6]],["mirror",[0,5,11]]]


    if Input_Crystal == "ORCC" or Input_Crystal == "oS":   
        [a,b,c,alpha,beta,gamma] = dim
        a1 = (a/2, -b/2, 0)
        a2 = (a/2, b/2, 0)
        a3 = (0, 0, c)
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        zeta = (1+a**2/b**2)/4
        Corners = np.zeros((10,3))
        Corners[0] = 0.0*b1+0.0*b2+0.0*b3
        Corners[1] = zeta*b1+zeta*b2+1/2*b3
        Corners[2] = -zeta*b1+(1-zeta)*b2+1/2*b3
        Corners[3] = 0.0*b1+1/2*b2+1/2*b3
        Corners[4] = 0.0*b1+1/2*b2+0.0*b3
        Corners[5] = -1/2*b1+1/2*b2+1/2*b3
        Corners[6] = zeta*b1+zeta*b2+0*b3
        Corners[7] = -zeta*b1+(1-zeta)*b2+0*b3
        Corners[8] = -1/2*b1+1/2*b2+0*b3
        Corners[9] = 0.0*b1+0.0*b2+1/2*b3
        Segments = np.array([[0,6],[0,8],[0,9],[1,3],[1,9],[2,5],[2,7],[3,4],[4,6],[5,8],[5,9],[7,8]])
        Plot_List = [0,6,4,3,1,9,0,8,7,2,5,8,"break",9,5]
        Corner_Labels = ['$\Gamma$','A','$A_1$','R','S','T','X','$X_1$','Y','Z']
        Ghosted_Segments = [[1,6],[2,3],[4,7]]
        Outer_Segments = [[2,5],[2,3],[3,1],[1,6],[7,2]]
        unfold_instructions = [["mirror",[0,6,1]],["mirror",[0,6,4]],["mirror",[0,8,5]]]
        

    if Input_Crystal == "HEX" or Input_Crystal == "hP":
        [a,b,c,alpha,beta,gamma] = dim
        a1 = (a/2, -(a*3**0.5)/2, 0)
        a2 = (a/2, (a*3**0.5)/2, 0)
        a3 = (0, 0, c)
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        Corners = np.zeros((6,3))
        Corners[0] = 0.0*b1+0.0*b2+0.0*b3
        Corners[1] = 0.0*b1+0.0*b2+1/2*b3
        Corners[2] = 1/3*b1+1/3*b2+1/2*b3
        Corners[3] = 1/3*b1+1/3*b2+0*b3
        Corners[4] = 1/2*b1+0.0*b2+1/2*b3
        Corners[5] = 1/2*b1+0.0*b2+0*b3
        Segments = np.array([[0,1],[0,3],[0,5],[1,2],[1,4],[2,3],[2,4],[3,5],[4,5]]) 
        Plot_List = [0,5,3,0,1,4,2,1,"break",4,5,"break",3,2]
        Corner_Labels = ['$\Gamma$','A','H','K','L','M']
        Ghosted_Segments = []
        Outer_Segments = [[4,2],[2,3]]
        unfold_instructions = [["mirror",[0,4,5]],["mirror",[0,3,5]],["mirror",[0,2,3]],\
                               ["2fmirror",[[0,5,4],[0,2,3]]],["2fmirror",[[0,2,3],[0,5,4]]]]

    
    if Input_Crystal=="RHL" or Input_Crystal=="hR":      
        [a,b,c,alpha,beta,gamma] = dim
        a1 = (a*np.cos(alpha/2), -a*np.sin(alpha/2), 0)
        a2 = (a*np.cos(alpha/2), a*np.sin(alpha/2), 0)
        a3 = (a*np.cos(alpha)/np.cos(alpha/2), 0, a*np.sqrt(1-(np.cos(alpha)**2)/(np.cos(alpha/2)**2)))     
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))

        if alpha < 90*np.pi/180:    #RHL1 
            eta = (1 + 4*np.cos(alpha))/(2 + 4*np.cos(alpha))
            nu = 3/4 - eta/2
            Corners = np.zeros((12,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = eta*b1+1/2*b2+(1-eta)*b3
            Corners[2] = 1/2*b1+(1-eta)*b2+(eta-1)*b3
            Corners[3] = 1/2*b1+1/2*b2+0*b3
            Corners[4] = 1/2*b1+0*b2+0*b3
            Corners[5] = 0.0*b1+0*b2+-1/2*b3
            Corners[6] = eta*b1+nu*b2+nu*b3
            Corners[7] = (1-nu)*b1+(1-nu)*b2+(1-eta)*b3
            Corners[8] = nu*b1+nu*b2+(eta-1)*b3
            Corners[9] = (1-nu)*b1+nu*b2+0*b3
            Corners[10] = nu*b1+0*b2+-nu*b3
            Corners[11] = 1/2*b1+1/2*b2+1/2*b3
            Segments = np.array([[0,4],[0,10],[0,11],[1,11],[2,4],[3,7],[3,9],[4,6],[7,11]])
            Plot_List = [0,4,2,"break",1,11,0,10,"break",9,3,7,11,"break",4,6]
            Corner_Labels = ["$\Gamma$",'B','$B_1$','F','L','$L_1$','P','$P_1$','$P_2$','Q','X','Z']
            Ghosted_Segments = [[11,6],[6,1],[1,7],[1,9],[9,2],[2,10],[2,8],[8,3],[8,5],[5,0],[5,10],[10,4]]
            Outer_Segments = [[10,2],[2,8],[2,1],[1,6],[1,7]]
            unfold_instructions = [["rot180",[10,0]],["mirror",[0,8,3]],["mirror",[0,4,6]],["2fmirror",[[0,3,8],[0,4,6]]]]
            print("RHL1")

        else:       #RHL2      #unsolved geometry problems, errors for alpha > 2ish radians
            eta = 1/(2*np.tan(alpha/2)**2)
            nu = 3/4 - eta/2
            Corners = np.zeros((8,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = 1/2*b1+-1/2*b2+0*b3
            Corners[2] = 1/2*b1+0*b2+0*b3
            Corners[3] = (1-nu)*b1-nu*b2+(1-nu)*b3
            Corners[4] = nu*b1+(nu-1)*b2+(nu-1)*b3
            Corners[5] = nu*b1+nu*b2+nu*b3
            Corners[6] = (1-nu)*b1+-nu*b2+-nu*b3
            Corners[7] = 1/2*b1+-1/2*b2+1/2*b3
            # Corners[4] = 2*Corners[6]-Corners[4]   #Sort of corrects screwy geometry?
            Segments = np.array([[0,1],[0,3],[0,5],[1,4],[2,6],[2,7],[3,7],[4,6],[5,7]])
            Plot_List = [0,3,7,5,0,1,4,6,2,7]
            Corner_Labels = ["$\Gamma$",'F','L','P','$P_1$','Q','$Q_1$','Z']
            Ghosted_Segments = [[5,2],[7,6],[4,0],[1,3]]
            Outer_Segments = [[5,7],[7,3],[7,6],[6,4]]
            unfold_instructions = [["rot180",[1,0]],["mirror",[0,7,3]],["mirror",[0,4,6]],["2fmirror",[[0,3,7],[0,4,6]]]]


    if Input_Crystal == "MCL" or Input_Crystal == "mP":   
        [a,b,c,alpha,beta,gamma] = dim
        a1 = (a, 0, 0)              #Ordering of a<=c, b<=c, alpha<pi/2, beta=gamma=pi/2
        a2 = (0, b, 0)
        a3 = (0, c*np.cos(alpha), c*np.sin(alpha))
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        eta = (1 - b*np.cos(alpha/c))/(2*np.sin(alpha)**2)
        nu = 1/2 - eta*c*np.cos(alpha/b)
        Corners = np.zeros((16,3))
        Corners[0] = 0.0*b1+0.0*b2+0.0*b3
        Corners[1] = 1/2*b1+1/2*b2+0*b3
        Corners[2] = 0*b1+1/2*b2+1/2*b3
        Corners[3] = 1/2*b1+0*b2+1/2*b3
        Corners[4] = 1/2*b1+0*b2+-1/2*b3
        Corners[5] = 1/2*b1+1/2*b2+1/2*b3
        Corners[6] = 0*b1+eta*b2+(1-nu)*b3
        Corners[7] = 0*b1+(1-eta)*b2+nu*b3
        Corners[8] = 0*b1+eta*b2+-nu*b3
        Corners[9] = 1/2*b1+eta*b2+(1-nu)*b3
        Corners[10] = 1/2*b1+(1-eta)*b2+nu*b3
        Corners[11] = 1/2*b1+eta*b2+-nu*b3
        Corners[12] = 0.0*b1+1/2*b2+0*b3
        Corners[13] = 0.0*b1+0.0*b2+1/2*b3
        Corners[14] = 0.0*b1+0.0*b2+-1/2*b3
        Corners[15] = 1/2*b1+0.0*b2+0*b3
        Segments = np.array([[0,13],[1,10],[1,12],[2,5],[2,6],[3,9],[3,13],[3,15],[5,10],[6,13],[7,12]])
        Plot_List = [0,13,6,2,5,10,1,12,7,"break",9,3,15,"break",13,3]
        Corner_Labels = ['$\Gamma$','A','C','D','$D_1$','E','H','$H_1$','$H_2$','M','$M_1$','$M_2$','X','Y','$Y_1$','Z']
        Ghosted_Segments = [[15,4],[4,14],[4,11],[11,8],[8,14],[8,12],[14,0],[11,1],[10,7],[7,2],[5,9],[9,6]]
        Outer_Segments = [[4,11],[11,8],[11,1],[1,10],[10,7],[10,5],[5,9],[9,6],[9,3]]
        unfold_instructions = [["mirror",[0,7,2]],["rot180",[0,15]]]


    if Input_Crystal=="MCLC" or Input_Crystal=="mS":      
        [a,b,c,alpha,beta,gamma] = dim      #Ordering of a<=c, b<=c, beta=gamma=pi/2
        a1 = (a/2, b/2, 0)
        a2 = (-a/2, b/2, 0)
        a3 = (0, c*np.cos(alpha), c*np.sin(alpha))
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        kgamma = np.arccos(np.dot(b1,b2)/(np.linalg.norm(b1)*np.linalg.norm(b2)))
        otherdecidingvar = b*np.cos(alpha/c) + b**2*np.sin(alpha/a**2)**2

        if modnumber == "2":    #MCLC2
            zeta = (2 - b*np.cos(alpha/c))/(4*np.sin(alpha)**2)
            eta = 2*zeta*c*np.cos(alpha/b)
            psi = 3/4 - a**2/(4*b**2*np.sin(alpha)**2)
            phi = psi + (3/4 - psi)*b*np.cos(alpha/c)
            Corners = np.zeros((15,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = 1/2*b1+0*b2+0*b3
            Corners[2] = 0*b1+-1/2*b2+0*b3
            Corners[3] = (1-zeta)*b1+(1-zeta)*b2+(1-eta)*b3    
            Corners[4] = zeta*b1+zeta*b2+eta*b3                   
            Corners[5] = -zeta*b1+-zeta*b2+(1-eta)*b3
            Corners[6] = (1-zeta)*b1+-zeta*b2+(1-eta)*b3
            Corners[7] = phi*b1+(1-phi)*b2+1/2*b3
            Corners[8] = (1-phi)*b1+(phi-1)*b2+1/2*b3
            Corners[9] = 1/2*b1+1/2*b2+1/2*b3
            Corners[10] = 1/2*b1+0*b2+1/2*b3
            Corners[11] = (1-psi)*b1+(psi-1)*b2+0*b3
            Corners[12] = 1/2*b1+1/2*b2+0*b3
            Corners[13] = -1/2*b1+-1/2*b2+0*b3
            Corners[14] = 0*b1+0*b2+1/2*b3          
            Segments = np.array([[0,1],[0,10],[0,12],[3,9],[3,12],[4,14],[7,9],[8,14]])
            Plot_List = [0,12,3,9,7,"break",8,14,4,"break",1,0,10]
            Corner_Labels = ["$\Gamma$",'N','$N_1$','F','$F_1$','$F_2$','$F_3$','I','$I_1$','L','M','X','Y','$Y_1$','Z']
            Ghosted_Segments = [[14,5],[5,8],[8,4],[4,9],[4,7],[7,3],[7,6],[6,8],[6,11],[11,1],[1,12],[11,2],[2,13],[13,0],[13,5]]
            Outer_Segments = [[12,3],[3,7],[7,4],[4,8],[8,5],[5,13],[7,6],[6,11],[6,8]]
            unfold_instructions = [["rot180",[11,0]],["mirror",[0,14,5]]]
            print("MLCL2222")
        
        elif kgamma > 90*np.pi/180:    #MCLC1 
            zeta = (2 - b*np.cos(alpha/c))/(4*np.sin(alpha)**2)
            eta = 2*zeta*c*np.cos(alpha/b)
            psi = 3/4 - a**2/(4*b**2*np.sin(alpha)**2)
            phi = psi + (3/4 - psi)*b*np.cos(alpha/c)
            Corners = np.zeros((16,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = 1/2*b1+0*b2+0*b3
            Corners[2] = 0*b1+-1/2*b2+0*b3
            Corners[3] = (1-zeta)*b1+(1-zeta)*b2+(1-eta)*b3    #tongue twister, because every code needs one of those
            Corners[4] = zeta*b1+zeta*b2+eta*b3                   #another
            Corners[5] = -zeta*b1+-zeta*b2+(1-eta)*b3
            Corners[6] = phi*b1+(1-phi)*b2+1/2*b3
            Corners[7] = (1-phi)*b1+(phi-1)*b2+1/2*b3
            Corners[8] = 1/2*b1+1/2*b2+1/2*b3
            Corners[9] = 1/2*b1+0*b2+1/2*b3
            Corners[10] = (1-psi)*b1+(psi-1)*b2+0*b3
            Corners[11] = psi*b1+(1-psi)*b2+0*b3
            Corners[12] = (psi-1)*b1+-psi*b2+0*b3
            Corners[13] = 1/2*b1+1/2*b2+0*b3
            Corners[14] = -1/2*b1+-1/2*b2+0*b3
            Corners[15] = 0*b1+0*b2+1/2*b3         
            Segments = np.array([[0,1],[0,9],[0,10],[0,13],[3,8],[3,13],[4,15],[6,8],[7,15],[11,13]])
            Plot_List = [0,13,3,8,6,"break",7,15,4,"break",13,11,"break",10,0,1,"break",9,0]
            Corner_Labels = ["$\Gamma$",'N','$N_1$','F','$F_1$','$F_2$','I','$I_1$','L','M','X','$X_1$','$X_2$','Y','$Y_1$','Z']
            Ghosted_Segments = []
            Outer_Segments = []
            unfold_instructions = []
            print("MLCL1")

        else:  
            if modnumber == "4":          #ORCF4    (Never actually triggers)  
                mu = (1+b**2/a**2)/4
                delta = b*c*np.cos(alpha/(2*a**2))
                zeta = mu - 1/4 + (1-b*np.cos(alpha/c))/(4*np.sin(alpha)**2)
                eta = 1/2 + 2*zeta*c*np.cos(alpha/b)
                phi = 1 + zeta - 2*mu
                psi = eta - 2*delta
                Corners = np.zeros((15,3))
                Corners[0] = 0.0*b1+0.0*b2+0.0*b3
                Corners[1] = (1-phi)*b1+(1-phi)*b2+(1-psi)*b3  
                Corners[2] = zeta*b1+zeta*b2+eta*b3       
                Corners[3] = (1-zeta)*b1+-zeta*b2+(1-eta)*b3  
                Corners[4] = -zeta*b1+-zeta*b2+(1-eta)*b3     
                Corners[5] = 1/2*b1+-1/2*b2+1/2*b3
                Corners[6] = 1/2*b1+0*b2+1/2*b3
                Corners[7] = 1/2*b1+0*b2+0*b3
                Corners[8] = 0*b1+-1/2*b2+0*b3
                Corners[9] = 1/2*b1+-1/2*b2+0*b3
                Corners[10] = mu*b1+mu*b2+delta*b3
                Corners[11] = (1-mu)*b1+-mu*b2+-delta*b3
                Corners[12] = -mu*b1+1-mu*b2+-delta*b3
                Corners[13] = mu*b1+(mu-1)*b2+delta*b3
                Corners[14] = 0*b1+0*b2+1/2*b3
                Segments = np.array([[0,6],[0,7],[0,9],[0,10],[1,2],[1,10],[2,14],[3,11],[5,14],[9,11]])
                Plot_List = [0,10,1,2,14,5,"break",3,11,9,0,7,"break",6,0]
                Corner_Labels = ["$\Gamma$",'F','H','$H_1$','$H_2$','I','M','N','$N_1$','X','Y','$Y_1$','$Y_2$','$Y_3$','Z']
                Ghosted_Segments = [[10,7],[7,11],[1,3],[3,5],[5,2],[5,4],[4,14],[4,12],[12,0],[12,8],[8,13],[13,5],[13,9]]
                Outer_Segments = [[10,1],[1,2],[2,5],[5,4],[4,12],[12,8],[8,13],[13,5],[5,3],[3,1],[3,11],[11,10]]
                unfold_instructions = [["rot180",[9,0]],["mirror",[0,14,4]]]
                print("MCLCfore")

            elif otherdecidingvar < 1:          #MCLC3 
                mu = (1+b**2/a**2)/4
                delta = b*c*np.cos(alpha/(2*a**2))
                zeta = mu - 1/4 + (1-b*np.cos(alpha/c))/(4*np.sin(alpha)**2)
                eta = 1/2 + 2*zeta*c*np.cos(alpha/b)
                phi = 1 + zeta - 2*mu
                psi = eta - 2*delta
                Corners = np.zeros((17,3))
                Corners[0] = 0.0*b1+0.0*b2+0.0*b3
                Corners[1] = (1-phi)*b1+(1-phi)*b2+(1-psi)*b3   #more tongue twisters!
                Corners[2] = phi*b1+(phi-1)*b2+psi*b3
                Corners[3] = (1-phi)*b1+-phi*b2+(1-psi)*b3  #not to be confused with 1-phi, 1-phi, 1-psi
                Corners[4] = zeta*b1+zeta*b2+eta*b3       
                Corners[5] = (1-zeta)*b1+-zeta*b2+(1-eta)*b3  
                Corners[6] = -zeta*b1+-zeta*b2+(1-eta)*b3      #guys, it's not funny anymore...
                Corners[7] = 1/2*b1+-1/2*b2+1/2*b3
                Corners[8] = 1/2*b1+0*b2+1/2*b3
                Corners[9] = 1/2*b1+0*b2+0*b3
                Corners[10] = 0*b1+-1/2*b2+0*b3
                Corners[11] = 1/2*b1+-1/2*b2+0*b3
                Corners[12] = mu*b1+mu*b2+delta*b3
                Corners[13] = (1-mu)*b1+-mu*b2+-delta*b3
                Corners[14] = -mu*b1+1-mu*b2+-delta*b3
                Corners[15] = mu*b1+(mu-1)*b2+delta*b3
                Corners[16] = 0*b1+0*b2+1/2*b3
                Segments = np.array([[0,8],[0,9],[0,11],[0,12],[1,4],[1,12],[2,7],[4,16],[5,13],[7,16],[11,13]])
                Plot_List = [0,12,1,4,16,7,2,"break",5,13,11,0,9,"break",8,0]
                Corner_Labels = ["$\Gamma$",'F','$F_1$','$F_2$','H','$H_1$','$H_2$','I','M','N','$N_1$','X','Y','Y$_1$','Y$_2$','Y$_3$','Z']
                Ghosted_Segments = [[16,6],[6,3],[3,7],[6,14],[3,15],[11,15],[15,14],[14,0],[4,2],[2,5],[5,1],[13,9],[9,12]]
                Outer_Segments = [[12,1],[1,4],[4,2],[2,7],[7,3],[3,6],[6,14],[15,14],[15,3],[2,5],[5,1],[5,13]]
                unfold_instructions = [["rot180",[11,0]],["mirror",[0,14,4]]]
                print("MCLC3 peoples")
                
            else:                   #MCLC5
                zeta = (b**2/a**2 + (1 - b*np.cos(alpha/c))/np.sin(alpha)**2)/4
                eta = 1/2 + 2*zeta*c*np.cos(alpha/b)
                rho = 1 - zeta*a**2/b**2
                mu = eta/2 + b**2/(4*a**2) - b*c*np.cos(alpha/(2*a**2))
                nu = 2*mu - zeta
                omega = (4*nu - 1 - b**2*np.sin(alpha/a**2)**2)*c/(2*b*np.cos(alpha))
                delta = zeta*c*np.cos(alpha/b) + omega/2 - 1/4
                Corners = np.zeros((19,3))
                Corners[0] = 0.0*b1+0.0*b2+0.0*b3
                Corners[1] = nu*b1+nu*b2+omega*b3  
                Corners[2] = (1-nu)*b1+(1-nu)*b2+(1-omega)*b3       
                Corners[3] = nu*b1+(nu-1)*b2+omega*b3    
                Corners[4] = zeta*b1+zeta*b2+eta*b3  
                Corners[5] = (1-zeta)*b1+-zeta*b2+(1-eta)*b3  
                Corners[6] = -zeta*b1+-zeta*b2+(1-eta)*b3  
                Corners[7] = rho*b1+(1-rho)*b2+1/2*b3  
                Corners[8] = (1-rho)*b1+(rho-1)*b2+1/2*b3     
                Corners[9] = 1/2*b1+-1/2*b2+1/2*b3
                Corners[10] = 1/2*b1+0*b2+1/2*b3
                Corners[11] = 1/2*b1+0*b2+0*b3
                Corners[12] = 0*b1+-1/2*b2+0*b3
                Corners[13] = 1/2*b1+-1/2*b2+0*b3
                Corners[14] = mu*b1+mu*b2+delta*b3
                Corners[15] = (1-mu)*b1+-mu*b2+-delta*b3
                Corners[16] = -mu*b1+1-mu*b2+-delta*b3
                Corners[17] = mu*b1+(mu-1)*b2+delta*b3
                Corners[18] = 0*b1+0*b2+1/2*b3
                Segments = np.array([[0,10],[0,11],[0,13],[0,14],[1,9],[1,14],[2,4],[4,18],[5,15],[7,9],[8,18],[13,15]])
                Plot_List = [0,14,1,9,7,"break",8,18,4,2,"break",5,15,13,0,11,"break",10,0]
                Corner_Labels = ["$\Gamma$",'F','$F_1$','$F_2$','H','$H_1$','$H_2$',\
                                 'I','$I_1$','L','M','N','$N_1$','X','Y','$Y_1$','$Y_2$','$Y_3$','Z']
                Ghosted_Segments = [[15,11],[11,14],[9,2],[2,7],[7,1],[7,5],[5,3],[3,8],[3,17],[17,13],[17,12],[12,16],[16,0],[16,6],[6,18],[6,8],[8,4]]
                Outer_Segments = [[14,1],[1,7],[7,2],[2,4],[4,8],[8,6],[6,16],[16,12],[12,17],[17,3],[3,8],[3,5],[5,7],[5,15]]
                unfold_instructions = [["rot180",[13,0]],["mirror",[0,14,4]]]
                print("MLCL5lves")


    if Input_Crystal=="TRI" or Input_Crystal=="aP":      
        [a,b,c,alpha,beta,gamma] = dim      
        a1 = (a, 0, 0)
        a2 = (b*np.cos(gamma), b*np.sin(gamma), 0)
        a3z = c/np.sin(gamma)*np.sqrt(np.sin(gamma)**2 - np.cos(alpha)**2 - np.cos(beta)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        a3 = (c*np.cos(beta), c/np.sin(gamma)*(np.cos(alpha)-np.cos(beta)*np.cos(gamma)), a3z)
        b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1,np.cross(a2,a3))
        b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2,np.cross(a3,a1))
        b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3,np.cross(a1,a2))
        kalpha = np.arccos(np.dot(b2,b3)/(np.linalg.norm(b2)*np.linalg.norm(b3)))
        
        if kalpha > 90*np.pi/180:    #TRI1a, TRI2a (TRI1a: kalpha, kbeta, kgamma > 90°) (TRI2a: kalpha, kbeta > 90°, kgamma = 90°)
            Corners = np.zeros((8,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = 1/2*b1+1/2*b2+0*b3
            Corners[2] = 0*b1+1/2*b2+1/2*b3 
            Corners[3] = 1/2*b1+0*b2+1/2*b3   
            Corners[4] = 1/2*b1+1/2*b2+1/2*b3                  
            Corners[5] = 1/2*b1+0*b2+0*b3
            Corners[6] = 0*b1+1/2*b2+0*b3
            Corners[7] = 0*b1+0*b2+1/2*b3      
            Segments = np.array([[0,1],[0,2],[0,3],[0,4],[0,5],[0,6],[0,7]])
            Plot_List = [5,0,6,"break",1,0,7,"break",3,0,2,"break",4,0]
            Corner_Labels = ["$\Gamma$",'L','M','N','R','X','Y','Z']
            Ghosted_Segments = []
            Outer_Segments = []
            unfold_instructions = []
            print("TRI1")

        else:    #TRI1b, TRI2b (TRI1b: kalpha, kbeta, kgamma < 90°) (TRI2b: kalpha, kbeta < 90°, kgamma = 90°)
            Corners = np.zeros((8,3))
            Corners[0] = 0.0*b1+0.0*b2+0.0*b3
            Corners[1] = 1/2*b1+-1/2*b2+0*b3
            Corners[2] = 0*b1+0*b2+1/2*b3 
            Corners[3] = -1/2*b1+-1/2*b2+1/2*b3   
            Corners[4] = 0*b1+-1/2*b2+1/2*b3                  
            Corners[5] = 0*b1+-1/2*b2+0*b3
            Corners[6] = 1/2*b1+0*b2+0*b3
            Corners[7] = -1/2*b1+0*b2+1/2*b3      
            Segments = np.array([[0,1],[0,2],[0,3],[0,4],[0,5],[0,6],[0,7]])
            Plot_List = [5,0,6,"break",1,0,7,"break",3,0,2,"break",4,0]
            Corner_Labels = ["$\Gamma$",'L','M','N','R','X','Y','Z']
            Ghosted_Segments = []
            Outer_Segments = []
            unfold_instructions = []
            print("TRI2222")



    #test for successful crystal conversion, and return to main program
    if b1[0] == None:
        print("Failure to recognize crystal ID code")
        exit()

    return (Corners,Segments,Plot_List,Corner_Labels,b1,b2,b3,Ghosted_Segments,Outer_Segments,unfold_instructions,a1,a2,a3)