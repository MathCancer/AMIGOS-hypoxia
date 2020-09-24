import numpy as np
import matplotlib.pyplot as plt


import cv2

######################################################################################################################  
# Test fraction of green cells for bias > 0.2 
def Read_FractionGreenCells(perct,initial_index = 0, last_index = 130):
    from pyMCDS import pyMCDS
    NumberGreen = np.zeros( (last_index-initial_index)+1 );
    NumberGreenResp = np.zeros( (last_index-initial_index)+1 );
    Fraction = np.zeros( (last_index-initial_index)+1 );
    for n in range( initial_index,last_index+1 ):
        filename='output'+"%08i"%n+'.xml'
        mcds=pyMCDS(filename,'output')
        geneY = mcds.data['discrete_cells']['genes_y']
        BiasMig = mcds.data['discrete_cells']['migration_bias'];
        cycle = mcds.data['discrete_cells']['cycle_model']
        biasDef = np.argwhere( (geneY > 0.5) & (BiasMig <= 0.2) & (cycle < 100) ).flatten()
        biasRes = np.argwhere( (geneY > 0.5) & (BiasMig > 0.2) & (cycle < 100) ).flatten()
        NumberGreen[n] = len(biasDef)+len(biasRes)
        NumberGreenResp[n] = len(biasRes)
        if (NumberGreen[n] != 0.0 ): Fraction[n] = NumberGreenResp[n]/NumberGreen[n]
    print(Fraction[-1])
    X = np.linspace(initial_index,last_index,(last_index-initial_index)+1)
    Ref = perct*np.ones(Fraction.shape[0])
    plt.plot(X,Fraction)
    plt.plot(X,Ref)
    plt.show()
    return

######################################################################################################################    
def ReadPositionGreenCells2Points():
    from astroML.correlation import two_point
    from astroML.correlation import bootstrap_two_point
    from pyMCDS import pyMCDS
    initial_index = 72;
    last_index = 72;
    for n in range( initial_index,last_index+1 ):
        filename='output'+"%08i"%n+'.xml'
        mcds=pyMCDS(filename,'output1')
        geneY = mcds.data['discrete_cells']['genes_y']
        cx = mcds.data['discrete_cells']['position_x'];
        cy = mcds.data['discrete_cells']['position_y'];
        cycle = mcds.data['discrete_cells']['cycle_model']
        GreenCells = np.argwhere( (geneY > 0.5) & (cx >= 0) & (cy >= 0) & (cycle < 100)).flatten()
        Position = np.array([cx[GreenCells],cy[GreenCells]]).T
        print(Position.shape)
        Twopointfunc(Position)
    return

def Twopointfunc(Position,NumBins=10):
    bins = np.linspace(0,500 , NumBins)
    corr, dcorr = bootstrap_two_point(Position, bins, Nbootstrap=5)
    print(corr)
    plt.plot(np.linspace(0,500 , NumBins-1),corr)

    plt.show()

######################################################################################################################  
# Read figure and classify according to Plumes, Escape cells, and Necrotic core 
def ImageOpenCV(file,plot=True):
    img = cv2.imread(file)
    blur = cv2.blur(img,(5,5))
    img = img[100:-70, 0:-1]
    img0 = img.copy()
    # Red filter
    channelred = img[:,:,0]
    _,thresh_val1 = cv2.threshold(channelred,15 ,255,cv2.THRESH_BINARY)
    thresh_val1 = cv2.bitwise_not(thresh_val1)
    # Green filter
    channelgreen = img[:,:,1]
    _,thresh_val2 = cv2.threshold(channelgreen,70 ,255,cv2.THRESH_BINARY)
    thresh_val2 = cv2.bitwise_not(thresh_val2)
    # Blue filter
    channelblue = img[:,:,2]
    _,thresh_val3 = cv2.threshold(channelblue,40 ,255,cv2.THRESH_BINARY)
    thresh_val3 = cv2.bitwise_not(thresh_val3)  
    # Combine
    mask = cv2.bitwise_and(thresh_val1,thresh_val2)
    mask2 = cv2.bitwise_and(thresh_val1,~mask)

    element = cv2.getStructuringElement(shape=cv2.MORPH_OPEN, ksize=(5, 5))
    element2 = cv2.getStructuringElement(shape=cv2.MORPH_OPEN, ksize=(8, 8))

    morph_img = thresh_val1.copy()
    cv2.morphologyEx(src=thresh_val1, op=cv2.MORPH_CLOSE, kernel=element, dst=morph_img)

    morph_img2 = mask2.copy()
    cv2.morphologyEx(src=mask2, op=cv2.MORPH_CLOSE, kernel=element, dst=morph_img2)

    morph_img3 = thresh_val3.copy()
    cv2.morphologyEx(src=thresh_val3, op=cv2.MORPH_CLOSE, kernel=element2, dst=morph_img3)

    contours,_ = cv2.findContours(morph_img,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    contours2,_ = cv2.findContours(morph_img2,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    contours3,_ = cv2.findContours(morph_img3,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)

    areas = [cv2.contourArea(c) for c in contours]
    sorted_areas = np.sort(areas)
    areas2 = [cv2.contourArea(c) for c in contours2]
    sorted_areas2 = np.sort(areas2)
    areas3 = [cv2.contourArea(c) for c in contours3]

    # cv2.imshow("RedMask", morph_img)
    # cv2.imshow("GreenMask", morph_img2)
    # cv2.imshow("BlueMask", morph_img3)
    # cv2.imwrite("onePic1.jpg", morph_img)
    # cv2.imwrite("onePic2.jpg", morph_img2)
    # cv2.imwrite("onePic3.jpg", morph_img3)
    # cv2.waitKey()
    # return
    
    
    #Ellipses
    img1 = Convert_BW_to_RGB(morph_img2)
    if (len(sorted_areas) > 0): 
        cntA=contours[areas.index(sorted_areas[-1])] #the biggest contour
        ellipseA = cv2.fitEllipse(cntA)
        cv2.ellipse(img,ellipseA,(0,0,0),2)    
    if (len(sorted_areas2) > 0): 
        cntA2=contours2[areas2.index(sorted_areas2[-1])] #the biggest contour
        epsilon = 0.00001 * cv2.arcLength(cntA2, True)
        approx = cv2.approxPolyDP(cntA2, epsilon, True)
        cv2.drawContours(img1, [approx], -1, (0, 255, 0), 4)
        ellipseA2 = cv2.fitEllipse(cntA2)
        cv2.ellipse(img1,ellipseA2,(100,100,0),2) 

    #Necrotic cells
    necrotic = False
    if( (np.sum(areas3)/np.sum(areas)) > 0.06):
        necrotic = True

    #Scapping cells
    scape = False  
    morph_img4 = morph_img2.copy()
    (h,k),(ma,MA),angle = ellipseA
    ma = ma + 0.05*ma
    MA = MA + 0.05*MA 
    ellipseB = (h,k),(ma,MA),angle 
    cv2.ellipse(morph_img4,ellipseB,(0,0,0),-1)
    img2 = Convert_BW_to_RGB(morph_img4)
    contours4,_ = cv2.findContours(morph_img4,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    areas4 = [cv2.contourArea(c) for c in contours4]
    if (np.sum(areas4) > 20.0):
        scape = True
    for ind in range(0,len(contours4)): 
        x = np.mean(contours4[ind][:,0,0])
        y = np.mean(contours4[ind][:,0,1])
        cv2.circle(img2,(int(x),int(y)),4,(0,200,0),2)

    #Plumes test
    plumes = False
    count = 0
    if (len(sorted_areas2) > 0): 
        (h,k),(ma,MA),angle = ellipseA2
        angle = angle*np.pi/180.0
    for i in range( 0,approx.shape[0]):
        (x,y) = approx[i,0,:]
        p = (((x - h)*np.cos(angle) + (y-k)*np.sin(angle))**2 / ((0.5*ma)**2)) +  (((x - h)*np.sin(angle) - (y-k)*np.cos(angle))**2 / ((0.5*MA)**2))
        if (p > 1.0): 
            dist = estimate_distance(x, y, 0.5*MA, 0.5*ma, x0=h, y0=k, angle=0, error=1e-5)
            tol =30.0
            if (dist >  tol):
                #print("Point: "+str(x)+" "+str(y)+" P: "+str(p)+" Dist: "+str(dist))
                #cv2.circle(img,(int(x),int(y)),4,(0,200,0),-1)
                cv2.circle(img1,(int(x),int(y)),4,(0,0,0),-1)
                plumes = True
    # t = np.linspace(0, 2*np.pi, 100)
    # plt.plot( h+0.5*ma*np.cos(t) , k+0.5*MA*np.sin(t) )
    # plt.plot(approx[:,0,0],approx[:,0,1])
    # plt.show()
    # if (plumes): print("There is plumes!")
    # else: print("There isn't plumes!")

    #Area Calculation
    # (h,k),(ma,MA),angle = ellipseA 
    # Area = (0.5*ma)*(0.5*MA)*np.pi
    # Fraction = 250.0/85.0
    # Real_Area = (Fraction**2)*Area
    # # Compatibilization
    # (x,y),radius = cv2.minEnclosingCircle(cntA)
    # print(radius)
    # center = (int(x),int(y))
    # radius = int(radius)
    # cv2.circle(img,center,radius,(0,0,0),2)     
    # print("Real ellipse area: "+str(Real_Area)+" micrometers")  

    if (plot): 
        #cv2.imshow("morph_img",morph_img)
        cv2.imshow("Plumes",img1)
        cv2.imshow("Escape cells",img2)
        cv2.imshow("Necrotic core",morph_img3)
        cv2.imshow("Original", img0)
        #cv2.imshow("img", img)
        #cv2.imshow("Green", GreenImage)
        # cv2.imwrite("Image.jpg", img0)
        # cv2.imwrite("Scaping.jpg", img)
        # cv2.imwrite("Plumes.jpg", GreenImage)
        cv2.waitKey()

    return plumes, scape, necrotic

def Convert_BW_to_RGB(image):
    New_image = np.zeros((image.shape[0],image.shape[1],3))
    ind = np.argwhere(image==0)
    New_image[ind[:,0],ind[:,1],0] = 0
    New_image[ind[:,0],ind[:,1],1] = 0
    New_image[ind[:,0],ind[:,1],2] = 0
    ind_0 = np.argwhere(image!=0)
    New_image[ind_0[:,0],ind_0[:,1],0] = 255
    New_image[ind_0[:,0],ind_0[:,1],1] = 255
    New_image[ind_0[:,0],ind_0[:,1],2] = 255
    return New_image
    
def Classify_all_parameter():
    Plumes = np.zeros(27)
    EscCell = np.zeros(27)
    NecCore = np.zeros(27)

    Bias = ["00","05","10"]
    Fraction = ["010","050","100"]
    TimePers = ["000","050","999"]

    index1 = []
    index2 = []
    index3 = []
    index4 = []
    index5 = []
    index6 = []
    index7 = []
    index8 = []
    index9 = []

    for n in range(len(Plumes)):
        if (n < 9):
            File = "output/Output_B"+Bias[0]
            index1.append(int(n))
        if ( (n>=9) & (n<18) ):
            File = "output/Output_B"+Bias[1]
            index2.append(int(n))
        if (n>=18):
            File = "output/Output_B"+Bias[2]
            index3.append(int(n))
        if (n%9 < 3):
            File += "_F"+Fraction[0]
            index4.append(int(n))
        if ( (n%9 >= 3) & (n%9 < 6) ):
            File += "_F"+Fraction[1]
            index5.append(int(n))
        if (n%9 >= 6):
            File += "_F"+Fraction[2]
            index6.append(int(n))
        if (n%3 == 0):
            File += "_T"+TimePers[0]
            index7.append(int(n))
        if (n%3 == 1):
            File += "_T"+TimePers[1]
            index8.append(int(n))
        if (n%3 == 2):
            File += "_T"+TimePers[2]
            index9.append(int(n))
        File += ".jpg"
        Plumes[n], EscCell[n], NecCore[n] = ImageOpenCV(File,False)
        print(str(n)+" - "+File +" "+str(Plumes[n])+" "+str(EscCell[n])+" "+str(NecCore[n]))
    # print(index1)
    # print(index2)
    # print(index3)
    # print(index4)
    # print(index5)
    # print(index6)
    # print(index7)
    # print(index8)
    # print(index9)

    # print("Bias: "+ Bias[0] + " Result"+ "-- Plumes:" +str(np.sum(Plumes[index1])/9.0) + " Escape: "+ str(np.sum(EscCell[index1])/9.0)+ " NecCore: "+ str(np.sum(NecCore[index1])/9.0))
    # print("Bias: "+ Bias[1] + " Result"+ "-- Plumes:" +str(np.sum(Plumes[index2])/9.0) + " Escape: "+ str(np.sum(EscCell[index2])/9.0)+ " NecCore: "+ str(np.sum(NecCore[index2])/9.0))
    # print("Bias: "+ Bias[2] + " Result"+ "-- Plumes:" +str(np.sum(Plumes[index3])/9.0) + " Escape: "+ str(np.sum(EscCell[index3])/9.0)+ " NecCore: "+ str(np.sum(NecCore[index3])/9.0))
    # print("Frac: "+ Fraction[0] + " Result"+ "-- Plumes:" +str(np.sum(Plumes[index4])/9.0) + " Escape: "+ str(np.sum(EscCell[index4])/9.0)+ " NecCore: "+ str(np.sum(NecCore[index4])/9.0))
    # print("Frac: "+ Fraction[1] + " Result"+ "-- Plumes:" +str(np.sum(Plumes[index5])/9.0) + " Escape: "+ str(np.sum(EscCell[index5])/9.0)+ " NecCore: "+ str(np.sum(NecCore[index5])/9.0))
    # print("Frac: "+ Fraction[2] + " Result"+ "-- Plumes:" +str(np.sum(Plumes[index6])/9.0) + " Escape: "+ str(np.sum(EscCell[index6])/9.0)+ " NecCore: "+ str(np.sum(NecCore[index6])/9.0))
    # print("Time: "+ TimePers[0] + " Result"+ "-- Plumes:" +str(np.sum(Plumes[index7])/9.0) + " Escape: "+ str(np.sum(EscCell[index7])/9.0)+ " NecCore: "+ str(np.sum(NecCore[index7])/9.0))
    # print("Time: "+ TimePers[1] + " Result"+ "-- Plumes:" +str(np.sum(Plumes[index8])/9.0) + " Escape: "+ str(np.sum(EscCell[index8])/9.0)+ " NecCore: "+ str(np.sum(NecCore[index8])/9.0))
    # print("Time: "+ TimePers[2] + " Result"+ "-- Plumes:" +str(np.sum(Plumes[index9])/9.0) + " Escape: "+ str(np.sum(EscCell[index9])/9.0)+ " NecCore: "+ str(np.sum(NecCore[index9])/9.0))  

    #generate data
    a = np.zeros((3, 9))
    a[0,0::3] = Plumes[index4[-3:]]
    a[0,1::3] = EscCell[index4[-3:]]
    a[0,2::3] = NecCore[index4[-3:]]

    a[1,0::3] = Plumes[index4[3:6]]
    a[1,1::3] = EscCell[index4[3:6]]
    a[1,2::3] = NecCore[index4[3:6]]

    a[2,0::3] = Plumes[index4[0:3]]
    a[2,1::3] = EscCell[index4[0:3]]
    a[2,2::3] = NecCore[index4[0:3]]
    discrete_matshow(a)

    b = np.zeros((3, 9))
    b[0,0::3] = Plumes[index5[-3:]]
    b[0,1::3] = EscCell[index5[-3:]]
    b[0,2::3] = NecCore[index5[-3:]]


    b[1,0::3] = Plumes[index5[3:6]]
    b[1,1::3] = EscCell[index5[3:6]]
    b[1,2::3] = NecCore[index5[3:6]]


    b[2,0::3] = Plumes[index5[0:3]]
    b[2,1::3] = EscCell[index5[0:3]]
    b[2,2::3] = NecCore[index5[0:3]]
    discrete_matshow(b)

    c = np.zeros((3, 9))
    c[0,0::3] = Plumes[index6[-3:]]
    c[0,1::3] = EscCell[index6[-3:]]
    c[0,2::3] = NecCore[index6[-3:]]

    c[1,0::3] = Plumes[index6[3:6]]
    c[1,1::3] = EscCell[index6[3:6]]
    c[1,2::3] = NecCore[index6[3:6]]

    c[2,0::3] = Plumes[index6[0:3]]
    c[2,1::3] = EscCell[index6[0:3]]
    c[2,2::3] = NecCore[index6[0:3]]

    discrete_matshow(c)
  
def discrete_matshow(data):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.colors import LinearSegmentedColormap
    fig, ax = plt.subplots(figsize=(6.5 ,7))
    #get discrete colormap
    #cmap = plt.get_cmap('RdGy', np.max(data)-np.min(data)+1)
    cm = LinearSegmentedColormap.from_list('MyBarColor', [(1.0,0,0),(0,0,1.0)], N=2)
    cmap = plt.get_cmap('RdGy', np.max(data)-np.min(data)+1)
    # set limits .5 outside true range
    mat = ax.matshow(data,cmap=cm,vmin = np.min(data)-.5, vmax = np.max(data)+.5,aspect='auto')
    #tell the colorbar to tick at integers
    divider = make_axes_locatable(ax)
    bar = divider.append_axes("right", size="5%", pad=0.1)
    cax = plt.colorbar(mat, ticks=[0,1], cax=bar)
    cax.ax.set_yticklabels(['False', 'True'],fontsize=14)
    ax.set_xticks([2.5,5.5])
    ax.set_xticklabels([])
    ax.set_yticks([0.5,1.5])
    ax.set_yticklabels([])
    ax.set_xticks([0.5,1.5,3.5,4.5,6.5,7.5], minor=True)
    plt.text(-39.8, 1.55, 'Plumes', fontsize=12,rotation=90)
    plt.text(-35.4, 1.55, 'Escape', fontsize=12,rotation=90)
    plt.text(-31, 1.55, 'Necrotic', fontsize=12,rotation=90)
    plt.text(-26.6, 1.55, 'Plumes', fontsize=12,rotation=90)
    plt.text(-22.2, 1.55, 'Escape', fontsize=12,rotation=90)
    plt.text(-17.8, 1.55, 'Necrotic', fontsize=12,rotation=90)
    plt.text(-13.4, 1.55, 'Plumes', fontsize=12,rotation=90)
    plt.text(-9.0, 1.55, 'Escape', fontsize=12,rotation=90)
    plt.text(-4.6, 1.55, 'Necrotic', fontsize=12,rotation=90)
    #T_p
    plt.text(-37.6, -0.58, '$T_p = 0h$', fontsize=14)
    plt.text(-24.95, -0.58, '$T_p = 50h$', fontsize=14)
    plt.text(-12.5, -0.58, '$T_p \geq 130h$', fontsize=14)
    #Bias \u03C6
    plt.text(-43.4, 1.05, '$\u03C6 = 1.0$', fontsize=14,rotation=90)
    plt.text(-43.4, 0.4, '$\u03C6 = 0.5$', fontsize=14,rotation=90)
    plt.text(-43.4, -0.25, '$\u03C6 = 0.1$', fontsize=14,rotation=90)
    ax.grid(which='minor', linewidth=1)
    ax.grid(color='w', linewidth=6)
    plt.show()  
  
###################################################################################################################### 

def Verify_symmetry(folder="outputStochastic"):
    from pyMCDS import pyMCDS
    index = 100
    subdiv = 10
    numfiles = 10
    filename = 'output'+"%08i"%index+'.xml'
    RedFrac = np.zeros((subdiv, numfiles))
    GreenFrac = np.zeros((subdiv, numfiles))
    NecFrac = np.zeros((subdiv, numfiles))
    for n in range(numfiles):
        path = folder+"/%02i"%n
        mcds=pyMCDS(filename,path)
        geneY = mcds.data['discrete_cells']['genes_y']
        cx = mcds.data['discrete_cells']['position_x'];
        cy = mcds.data['discrete_cells']['position_y'];
        cycle = mcds.data['discrete_cells']['cycle_model']
        GreenCells = np.argwhere( (geneY > 0.5) & (cycle < 100)).flatten()
        RedCells = np.argwhere( (geneY < 0.5) & (cycle < 100)).flatten()
        NecCells = np.argwhere((cycle > 100)).flatten()
        theta = np.arctan((cy/cx))
        theta[np.isnan(theta)] = 0.0
        quadrant1 = np.argwhere( (cx >= 0.0) & (cy >= 0.0)).flatten()
        quadrant2 = np.argwhere( (cx < 0.0) & (cy >= 0.0)).flatten()
        quadrant3 = np.argwhere( (cx <= 0.0) & (cy < 0.0)).flatten()
        quadrant4 = np.argwhere( (cx > 0.0) & (cy < 0.0)).flatten()
        theta[quadrant2] += np.pi
        theta[quadrant3] += np.pi
        theta[quadrant4] += 2*np.pi
        theta0 = 0.0
        thetaf = (2*np.pi/subdiv)
        for angle in range(subdiv):
            indR = np.argwhere( (theta >= theta0) & (theta < thetaf) & (geneY < 0.5) & (cycle < 100)).flatten()
            indG = np.argwhere( (theta >= theta0) & (theta < thetaf) & (geneY > 0.5) & (cycle < 100)).flatten()
            indN = np.argwhere( (theta >= theta0) & (theta < thetaf) & (cycle > 100)).flatten()
            Numcells = len(RedCells) + len(GreenCells) + len(NecCells)
            if (len(RedCells) > 0): RedFrac[angle,n] = len(indR)/Numcells
            if (len(GreenCells) > 0): GreenFrac[angle,n] = len(indG)/Numcells
            if (len(NecCells) > 0): NecFrac[angle,n] = len(indN)/Numcells
            # plt.figure()
            # plt.scatter(cx[indR],cy[indR],color='red')
            # plt.scatter(cx[indG],cy[indG],color='green')
            # plt.scatter(cx[indN],cy[indN],color='blue')
            # x = np.linspace(-1500.0,1500.0,100)
            # plt.plot(x,np.tan(theta0)*x)
            # plt.plot(x,np.tan(thetaf)*x)
            # plt.xlim(-1500.0,1500.0)
            # plt.ylim(-1500.0,1500.0)
            theta0 = thetaf
            thetaf = thetaf + (2*np.pi/subdiv)
        # plt.show()
    # Plot_Symmetry(np.mean(RedFrac),RedFrac[:,0],'#cc0000')
    # Plot_Symmetry(np.mean(GreenFrac),GreenFrac[:,0],'#007f00')
    # Plot_Symmetry(np.mean(NecFrac),NecFrac[:,3],'#6336de')
    return RedFrac,GreenFrac,NecFrac
  
def Plot_Symmetry(mean,Fraction,color):
    from matplotlib.patches import Circle, Wedge
    from matplotlib.collections import PatchCollection
    theta = np.linspace(0, 2*np.pi, 100)
    r = mean
    fig, ax = plt.subplots()
    origin = np.zeros(2)
    patches = []
    theta0 = 0
    thetaf = (360.0/Fraction.shape[0])
    for i in range( 0,Fraction.shape[0] ):
        wedge = Wedge(origin,  Fraction[i], theta0, thetaf)
        patches.append(wedge)
        theta0 = thetaf
        thetaf = thetaf + (360.0/Fraction.shape[0])
  
    #colors = np.linspace(0,len(patches),len(patches))
    p = PatchCollection(patches, alpha=0.8)
    p.set_color(color)
    ax.add_collection(p)
    #Circle
    c = Circle(origin, r, fill=False, edgecolor='gray',ls=':')
    p1 = PatchCollection([c], match_original=True)
    ax.add_collection(p1)
    plt.xlim(-1.0*mean - 0.5*mean,1.0*mean + 0.5*mean)
    plt.ylim(-1.0*mean - 0.5*mean,1.0*mean + 0.5*mean)
    ax.set_aspect(1)

    plt.show()

def Plot_Symmetry_Hist(Fraction1,Fraction2,Fraction3):
    origin = np.zeros(2)
    patches1 = []
    patches2 = []
    patches3 = []
    opacity = 0.5
  
    #Red
    theta0 = 0
    thetaf = (360.0/Fraction3.shape[0])
    for i in range( 0,Fraction1.shape[0] ):
        wedge = Wedge(origin,np.mean(Fraction1[i,:]), theta0, thetaf)
        patches1.append(wedge)
        theta0 = thetaf
        thetaf = thetaf + (360.0/Fraction1.shape[0])
    p1 = PatchCollection(patches1, alpha=opacity)
    p1.set_color('#cc0000')
  
    # Green
    theta0 = 0
    thetaf = (360.0/Fraction3.shape[0])
    for i in range( 0,Fraction2.shape[0] ):
        wedge = Wedge(origin, np.mean(Fraction2[i,:]), theta0, thetaf)
        patches2.append(wedge)
        theta0 = thetaf
        thetaf = thetaf + (360.0/Fraction2.shape[0])
    p2 = PatchCollection(patches2, alpha=opacity)
    p2.set_color('#007f00')
  
    # Necrotic
    theta0 = 0
    thetaf = (360.0/Fraction3.shape[0])
    for i in range( 0,Fraction3.shape[0] ):
        wedge = Wedge(origin, np.mean(Fraction3[i,:]), theta0, thetaf)
        patches3.append(wedge)
        theta0 = thetaf
        thetaf = thetaf + (360.0/Fraction3.shape[0])
    p3 = PatchCollection(patches3, alpha=opacity)
    p3.set_color('#6336de')
  
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(1,2,1)
    ax.set_title('Average fraction for each slice')
    tics = np.linspace(-0.10,0.10,num=5)
    ax.set_xticks(tics)
    ax.set_yticks(tics)
    ax.set_xlim(tics[0],tics[-1])
    ax.set_ylim(tics[0],tics[-1])
    ax.ticklabel_format(style='sci',axis='both')
    p = [p1,p2,p3]
    Average = np.array([np.mean(Fraction1),np.mean(Fraction2),np.mean(Fraction3)])
    ind = np.argsort(Average)
    ax.add_collection(p[ind[2]])
    ax.add_collection(p[ind[1]])
    ax.add_collection(p[ind[0]])
    # Label
    #LabelRadius = np.max(Average) + 0.01
    Fraction = [Fraction1,Fraction2,Fraction3]
    LabelRadius = np.mean(Fraction[ind[2]],axis=1) + 0.01
    LabelTheta = np.linspace((np.pi/Fraction3.shape[0]),2*np.pi-(np.pi/Fraction3.shape[0]),10)
    LabelPosition = np.zeros((Fraction3.shape[0],2))
    LabelPosition[:,0] = LabelRadius*np.cos(LabelTheta)
    LabelPosition[:,1] = LabelRadius*np.sin(LabelTheta)
    Label = []
    [Label.append("%i"%i) for i in range(1,Fraction3.shape[0] + 1)]
    for i in range(Fraction3.shape[0]):
        ax.text(LabelPosition[i,0] , LabelPosition[i,1], Label[i], fontsize=12, color='gray')
    ax.set_aspect(1)

    x = np.arange(1,Fraction1.shape[0]+1)
    ax2 = fig.add_subplot(3,2,2)
    ax2.errorbar(x,np.mean(Fraction1,axis=1),yerr=np.std(Fraction1,axis=1),fmt='o',color='#cc0000')#'#7a0000'
    # ax2.errorbar(x,np.mean(AFraction1,axis=1),yerr=np.std(Fraction1,axis=1),fmt='o',color='#a30000')
    # ax2.errorbar(x,np.mean(BFraction1,axis=1),yerr=np.std(Fraction1,axis=1),fmt='o',color='#cc0000')
    ax2.set_xticks(x)
    ax2.set_xlabel('Slices')
    ax2.set_ylabel('Fraction')
    ax2.ticklabel_format(style='sci',scilimits=(-2,4),axis='y',useMathText=True)

    ax3 = fig.add_subplot(3,2,4)
    ax3.errorbar(x,np.mean(Fraction2,axis=1),yerr=np.std(Fraction2,axis=1),fmt='o',color='#007f00')#'#001900'
    # ax3.errorbar(x,np.mean(AFraction2,axis=1),yerr=np.std(Fraction2,axis=1),fmt='o',color='#004c00')
    # ax3.errorbar(x,np.mean(BFraction2,axis=1),yerr=np.std(Fraction2,axis=1),fmt='o',color='#007f00')
    ax3.set_xticks(x)
    ax3.set_xlabel('Slices')
    ax3.set_ylabel('Fraction')
    ax3.ticklabel_format(style='sci',scilimits=(-2,4),axis='y',useMathText=True)

    ax4 = fig.add_subplot(3,2,6)
    ax4.errorbar(x,np.mean(Fraction3,axis=1),yerr=np.std(Fraction3,axis=1),fmt='o',color='#6336de')
    # ax4.errorbar(x,np.mean(AFraction3,axis=1),yerr=np.std(Fraction3,axis=1),fmt='o',color='#6336de')
    # ax4.errorbar(x,np.mean(BFraction3,axis=1),yerr=np.std(Fraction3,axis=1),fmt='o',color='#6336de')
    ax4.set_xticks(x)
    ax4.set_xlabel('Slices')
    ax4.set_ylabel('Fraction')
    ax4.ticklabel_format(style='sci',scilimits=(-2,4),axis='y',useMathText=True)

    plt.subplots_adjust(wspace=0.4,hspace=0.5)
    plt.show()  

def Run_all_cases():
    RedFrac1,GreenFrac1,NecFrac1 = Verify_symmetry("outputStochastic1")
    RedFrac2,GreenFrac2,NecFrac2 = Verify_symmetry("outputStochastic05")
    RedFrac3,GreenFrac3,NecFrac3 = Verify_symmetry("outputStochastic01791")
    Plot_Symmetry_Hist(RedFrac1,GreenFrac1,NecFrac1)
    Plot_Symmetry_Hist(RedFrac2,GreenFrac2,NecFrac2)
    Plot_Symmetry_Hist(RedFrac3,GreenFrac3,NecFrac3)
###################################################################################################################### 
 
def ReplicaAnalysis(folder="OneCase"):
    Plumes = np.zeros(20)
    EscCell = np.zeros(20)
    NecCore = np.zeros(20)
  
    for n in range(len(Plumes)):
        File = folder+"/final"+"%02i"%(n+1)+".jpg"
        Plumes[n], EscCell[n], NecCore[n] = ImageOpenCV(File,False)
        print(str(n+1)+" Plumes: " + str(Plumes[n]) + " Escape: " + str(EscCell[n]) + " NecCore: " + str(NecCore[n]))
    
    PlumesT = 100.0*np.sum(Plumes)/len(Plumes)
    PlumesF = 100.0 - PlumesT
    EscapeT = 100.0*np.sum(EscCell)/len(Plumes)
    EscapeF = 100.0 - EscapeT
    NecCorT = 100.0*np.sum(NecCore)/len(Plumes)
    NecCorF = 100.0 - NecCorT
    print(">>> Plumes: " + str(PlumesT) + "% Escape: " + str(EscapeT) + "% NecCore: " + str(NecCorT)+"%")
    
    labels = 'True', 'False'
    colors = ['Blue', 'Red']
    sizesP = [PlumesT, PlumesF]
    sizesE = [EscapeT, EscapeF]
    sizesN = [NecCorT, NecCorF]
    explode = (0, 0.5)  # only "explode" the 2nd slice (i.e. 'Hogs')
    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(1,3,1)
    ax1.pie(sizesP, explode=explode, autopct='%1.0f%%', shadow=True, startangle=90,colors=colors,textprops={'fontsize': 16})
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax2 = fig.add_subplot(1,3,2)
    ax2.pie(sizesE, explode=explode, autopct='%1.0f%%', shadow=True, startangle=90,colors=colors,textprops={'fontsize': 16})
    ax2.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax3 = fig.add_subplot(1,3,3)
    ax3.pie(sizesN, explode=explode, autopct='%1.0f%%', shadow=True, startangle=90,colors=colors,textprops={'fontsize': 16})
    ax3.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    fig.legend(labels, loc="right",prop={'size': 16})
    plt.show()
  
######################################################################################################################     
     
 
from math import sin, cos, atan2, pi, fabs
def ellipe_tan_dot(rx, ry, px, py, theta):
    '''Dot product of the equation of the line formed by the point
    with another point on the ellipse's boundary and the tangent of the ellipse
    at that point on the boundary.
    '''
    return ((rx ** 2 - ry ** 2) * cos(theta) * sin(theta) -
            px * rx * sin(theta) + py * ry * cos(theta))


def ellipe_tan_dot_derivative(rx, ry, px, py, theta):
    '''The derivative of ellipe_tan_dot.
    '''
    return ((rx ** 2 - ry ** 2) * (cos(theta) ** 2 - sin(theta) ** 2) -
            px * rx * cos(theta) - py * ry * sin(theta))


def estimate_distance(x, y, rx, ry, x0=0, y0=0, angle=0, error=1e-5):
    '''Given a point (x, y), and an ellipse with major - minor axis (rx, ry),
    its center at (x0, y0), and with a counter clockwise rotation of
    `angle` degrees, will return the distance between the ellipse and the
    closest point on the ellipses boundary.
    '''
    x -= x0
    y -= y0
    if angle:
        # rotate the points onto an ellipse whose rx, and ry lay on the x, y
        # axis
        angle = -pi / 180. * angle
        x, y = x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle)

    theta = atan2(rx * y, ry * x)
    while fabs(ellipe_tan_dot(rx, ry, x, y, theta)) > error:
        theta -= ellipe_tan_dot(
            rx, ry, x, y, theta) / \
            ellipe_tan_dot_derivative(rx, ry, x, y, theta)

    px, py = rx * cos(theta), ry * sin(theta)
    return ((x - px) ** 2 + (y - py) ** 2) ** .5
