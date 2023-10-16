import xarray as xr
import numpy as np
def window3x3_a(arr, shape=(3, 3)):
    x,y = np.shape(arr)
    
    xmin_a=np.zeros(x).astype('i')
    xmax_a=np.zeros(x).astype('i')
    ymin_a=np.zeros(y).astype('i')
    ymax_a=np.zeros(y).astype('i')
    
    
    
   
    
    r_win = np.floor(shape[0] / 2).astype(int)
    c_win = np.floor(shape[1] / 2).astype(int)
    for i in range(x):
        xmin = max(0, i - r_win)
        xmax = min(x, i + r_win + 1)
        xmin_a[i]=xmin
        xmax_a[i]=xmax
    
        for j in range(y):
            ymin = max(0, j - c_win)
            ymax = min(y, j + c_win + 1)
            ymin_a[j]=ymin
            ymax_a[j]=ymax
       
    return(xmin_a,xmax_a,ymin_a,ymax_a)
if CREATE:
def create_bed_slope_angle(d4)
""" d4 = dimarry with bed vairable , 
    return: returns an array with the shape of the grid with the computed slopes
    I am assumin equally 1000m girdsize


    The algorithm calculates the gradient using a first-order forward or backward difference on the corner points, first
    order central differences at the boarder points, and a 3x3 moving window for every cell with 8 surrounding cells (in
    the middle of the grid) using a third-order finite difference weighted by reciprocal of squared distance

    Assumed 3x3 window:

                        -------------------------
                        |   a   |   b   |   c   |
                        -------------------------
                        |   d   |   e   |   f   |
                        -------------------------
                        |   g   |   h   |   i   |
                        -------------------------


    """"
    
    bed = d4['bed'].values.copy()
    mnan = np.isnan(bed)
    bed[mnan]= -9999.0
            
    x = d3.x.values.copy()
    y = d3.y.values.copy()
    nx= np.size(x)
    ny= np.size(y)
    dx = 1000 #made it 1km eesolution in matlab script
    dy=1000 #made it 1km eesolution in matlab script
    # gen = window3x3(bd)
    # windows_3x3 = np.asarray(list(gen))
    xmin_a,xmax_a,ymin_a,ymax_a = window3x3_a(bed)

    dzdx = np.empty((ny, nx))
    dzdy = np.empty((ny, nx))
    loc_string = np.empty((ny, nx), dtype="S25")

    for ax_y in range(ny):
            for ax_x in range(nx):
                windows = bed[xmin_a[ax_y]:xmax_a[ax_y],ymin_a[ax_x]:ymax_a[ax_x]]

                # corner points
                if ax_x == 0 and ax_y == 0:  # top left corner
                    dzdx[ax_y, ax_x] = (windows[0][1] - windows[0][0]) / dx

                    dzdy[ax_y, ax_x] = (windows[1][0] - windows[0][0]) / dy
                    loc_string[ax_y, ax_x] = 'top left corner'

                elif ax_x == nx - 1 and ax_y == 0:  # top right corner
                    dzdx[ax_y, ax_x] = (windows[0][1] - windows[0][0]) / dx
                    dzdy[ax_y, ax_x] = (windows[1][1] - windows[0][1]) / dy
                    loc_string[ax_y, ax_x] = 'top right corner'

                elif ax_x == 0 and ax_y == ny - 1:  # bottom left corner
                    dzdx[ax_y, ax_x] = (windows[1][1] - windows[1][0]) / dx
                    dzdy[ax_y, ax_x] = (windows[1][0] - windows[0][0]) / dy
                    loc_string[ax_y, ax_x] = 'bottom left corner'

                elif ax_x == nx - 1 and ax_y == ny - 1:  # bottom right corner
                    dzdx[ax_y, ax_x] = (windows[1][1] - windows[1][0]) / dx
                    dzdy[ax_y, ax_x] = (windows[1][1] - windows[0][1]) / dy
                    loc_string[ax_y, ax_x] = 'bottom right corner'

                # top boarder
                elif (ax_y == 0) and (ax_x != 0 and ax_x != nx - 1):
                    dzdx[ax_y, ax_x] = (windows[0][-1] - windows[0][0]) / (2 * dx)
                    dzdy[ax_y, ax_x] = (windows[1][1] - windows[0][1]) / dy
                    loc_string[ax_y, ax_x] = 'top boarder'

                # bottom boarder
                elif ax_y == ny - 1 and (ax_x != 0 and ax_x != nx - 1):
                    dzdx[ax_y, ax_x] = (windows[1][-1] - windows[1][0]) / (2 * dx)
                    dzdy[ax_y, ax_x] = (windows[1][1] - windows[0][1]) / dy
                    loc_string[ax_y, ax_x] = 'bottom boarder'

                # left boarder
                elif ax_x == 0 and (ax_y != 0 and ax_y != ny - 1):
                    dzdx[ax_y, ax_x] = (windows[1][1] - windows[1][0]) / dx
                    dzdy[ax_y, ax_x] = (windows[-1][0] - windows[0][0]) / (2 * dy)
                    loc_string[ax_y, ax_x] = 'left boarder'

                # right boarder
                elif ax_x == nx - 1 and (ax_y != 0 and ax_y != ny - 1):
                    dzdx[ax_y, ax_x] = (windows[1][1] - windows[1][0]) / dx
                    dzdy[ax_y, ax_x] = (windows[-1][-1] - windows[0][-1]) / (2 * dy)
                    loc_string[ax_y, ax_x] = 'right boarder'

                # middle grid
                else:
                    a = windows[0][0]
                    b = windows[0][1]
                    c = windows[0][-1]
                    d = windows[1][0]
                    f = windows[1][-1]
                    g = windows[-1][0]
                    h = windows[-1][1]
                    i = windows[-1][-1]

                    dzdx[ax_y, ax_x] = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * dx)
                    dzdy[ax_y, ax_x] = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * dy)
                    loc_string[ax_y, ax_x] = 'middle grid'
    hpot = np.hypot(abs(dzdy), abs(dzdx))
    slopes_angle = np.degrees(np.arctan(hpot))
    slopes_angle_ana =slopes_angle.copy()
    slopes_angle_ana[slopes_angle_ana >70]=np.nan
    ds = xr.Dataset(
    data_vars=dict(
        bed_slope=(["x", "y", ], slopes_angle_ana )))
    ds.to_netcdf('./../Data/Fields/bedslope.nc',mode='w')
    return(slopes_angle_ana)


def window3x3_a(arr, shape=(3, 3)):
    x,y = np.shape(arr)
    
    xmin_a=np.zeros(x).astype('i')
    xmax_a=np.zeros(x).astype('i')
    ymin_a=np.zeros(y).astype('i')
    ymax_a=np.zeros(y).astype('i')
    
    
    
   
    
    r_win = np.floor(shape[0] / 2).astype(int)
    c_win = np.floor(shape[1] / 2).astype(int)
    for i in range(x):
        xmin = max(0, i - r_win)
        xmax = min(x, i + r_win + 1)
        xmin_a[i]=xmin
        xmax_a[i]=xmax
    
        for j in range(y):
            ymin = max(0, j - c_win)
            ymax = min(y, j + c_win + 1)
            ymin_a[j]=ymin
            ymax_a[j]=ymax
       
    return(xmin_a,xmax_a,ymin_a,ymax_a)
            

def nnmean(mat):

    """

    :mat =matrix has to be ar leest 3 x 
   


    Assumed 3x3 window:

                        -------------------------
                        |   a   |   b   |   c   |
                        -------------------------
                        |   d   |   e   |   f   |
                        -------------------------
                        |   g   |   h   |   i   |
                        -------------------------


    """
#     x = d3.x.values.copy()
#     y = d3.y.values.copy()
    ny,nx = mat.shape
   
    # gen = window3x3(bd)
    # windows_3x3 = np.asarray(list(gen))
    xmin_a,xmax_a,ymin_a,ymax_a = window3x3_a(mat)

    dzdx = np.empty((ny, nx)) #matrix for nanmean
  
    loc_string = np.empty((ny, nx), dtype="S25")
    for ax_y in range(ny):
            for ax_x in range(nx):
                windows = mat[xmin_a[ax_y]:xmax_a[ax_y],ymin_a[ax_x]:ymax_a[ax_x]]

                # corner points
                if ax_x == 0 and ax_y == 0:  # top left corner
                   
                    dzdx[ax_y, ax_x] = np.nanmean([windows[0][1] ,windows[1][0] , windows[1][1]]) 

                    
                    loc_string[ax_y, ax_x] = 'top left corner'

                elif ax_x == nx - 1 and ax_y == 0:  # top right corner
                    dzdx[ax_y, ax_x] = np.nanmean([windows[1][1],windows[0][0],windows[1][0]])
               
                    loc_string[ax_y, ax_x] = 'top right corner'

                elif ax_x == 0 and ax_y == ny - 1:  # bottom left corner
                    dzdx[ax_y, ax_x] = np.nanmean([windows[1][1],windows[0][0],windows[0][1]])
                
                    loc_string[ax_y, ax_x] = 'bottom left corner'

                elif ax_x == nx - 1 and ax_y == ny - 1:  # bottom right corner
                    dzdx[ax_y, ax_x] = np.nanmean([windows[1][0],windows[0][0],windows[0][1]])
                    loc_string[ax_y, ax_x] = 'bottom right corner'

                # top boarder
                elif (ax_y == 0) and (ax_x != 0 and ax_x != nx - 1):
                    dzdx[ax_y, ax_x] = np.nanmean([windows[0][0],windows[0][-1],windows[1][0],windows[1][1],windows[1][2]])
                    loc_string[ax_y, ax_x] = 'top boarder'

                # bottom boarder
                elif ax_y == ny - 1 and (ax_x != 0 and ax_x != nx - 1):
                    dzdx[ax_y, ax_x] = np.nanmean([windows[1][0],windows[1][-1],windows[0][0],windows[0][1],windows[0][2]])
                    
                    loc_string[ax_y, ax_x] = 'bottom boarder'

                # left boarder
                elif ax_x == 0 and (ax_y != 0 and ax_y != ny - 1):
                    dzdx[ax_y, ax_x] = np.nanmean([windows[0][0],windows[1][-1],windows[2][0],windows[0][1],windows[2][1]])
                    loc_string[ax_y, ax_x] = 'left boarder'

                # right boarder
                elif ax_x == nx - 1 and (ax_y != 0 and ax_y != ny - 1):
                    dzdx[ax_y, ax_x] = np.nanmean([windows[0][0],windows[2][0],windows[0][1],windows[2][1],windows[1][0],])
                   
                    loc_string[ax_y, ax_x] = 'right boarder'

                # middle grid
                else:
                    a = windows[0][0]
                    b = windows[0][1]
                    c = windows[0][-1]
                    d = windows[1][0]
                    f = windows[1][-1]
                    g = windows[-1][0]
                    h = windows[-1][1]
                    i = windows[-1][-1]

                    dzdx[ax_y, ax_x] = np.nanmean([a,b,c,d,f,g,h,i])
                    
                    loc_string[ax_y, ax_x] = 'middle grid'
    return(dzdx,loc_string)



def nnmean_skipnan(mat):
    '''mat =matrix has to be ar least 3 x 3'''
#     x = d3.x.values.copy()
#     y = d3.y.values.copy()
    ny,nx = mat.shape
   
    # gen = window3x3(bd)
    # windows_3x3 = np.asarray(list(gen))
    xmin_a,xmax_a,ymin_a,ymax_a = window3x3_a(mat)

    dzdx = np.empty((ny, nx)) #matrix for nanmean
  
    loc_string = np.empty((ny, nx), dtype="S25")
    for ax_y in range(ny):
            for ax_x in range(nx):
                if np.isnan(mat[ax_y, ax_x]):
                    dzdx[ax_y, ax_x] =np.nan
                    continue
                
                windows = mat[xmin_a[ax_y]:xmax_a[ax_y],ymin_a[ax_x]:ymax_a[ax_x]]

                # corner points
                if ax_x == 0 and ax_y == 0:  # top left corner
                   
                    dzdx[ax_y, ax_x] = np.nanmean([windows[0][1] ,windows[1][0] , windows[1][1]]) 

                    
                    loc_string[ax_y, ax_x] = 'top left corner'

                elif ax_x == nx - 1 and ax_y == 0:  # top right corner
                    dzdx[ax_y, ax_x] = np.nanmean([windows[1][1],windows[0][0],windows[1][0]])
               
                    loc_string[ax_y, ax_x] = 'top right corner'

                elif ax_x == 0 and ax_y == ny - 1:  # bottom left corner
                    dzdx[ax_y, ax_x] = np.nanmean([windows[1][1],windows[0][0],windows[0][1]])
                
                    loc_string[ax_y, ax_x] = 'bottom left corner'

                elif ax_x == nx - 1 and ax_y == ny - 1:  # bottom right corner
                    dzdx[ax_y, ax_x] = np.nanmean([windows[1][0],windows[0][0],windows[0][1]])
                    loc_string[ax_y, ax_x] = 'bottom right corner'

                # top boarder
                elif (ax_y == 0) and (ax_x != 0 and ax_x != nx - 1):
                    dzdx[ax_y, ax_x] = np.nanmean([windows[0][0],windows[0][-1],windows[1][0],windows[1][1],windows[1][2]])
                    loc_string[ax_y, ax_x] = 'top boarder'

                # bottom boarder
                elif ax_y == ny - 1 and (ax_x != 0 and ax_x != nx - 1):
                    dzdx[ax_y, ax_x] = np.nanmean([windows[1][0],windows[1][-1],windows[0][0],windows[0][1],windows[0][2]])
                    
                    loc_string[ax_y, ax_x] = 'bottom boarder'

                # left boarder
                elif ax_x == 0 and (ax_y != 0 and ax_y != ny - 1):
                    dzdx[ax_y, ax_x] = np.nanmean([windows[0][0],windows[1][-1],windows[2][0],windows[0][1],windows[2][1]])
                    loc_string[ax_y, ax_x] = 'left boarder'

                # right boarder
                elif ax_x == nx - 1 and (ax_y != 0 and ax_y != ny - 1):
                    dzdx[ax_y, ax_x] = np.nanmean([windows[0][0],windows[2][0],windows[0][1],windows[2][1],windows[1][0],])
                   
                    loc_string[ax_y, ax_x] = 'right boarder'

                # middle grid
                else:
                    a = windows[0][0]
                    b = windows[0][1]
                    c = windows[0][-1]
                    d = windows[1][0]
                    f = windows[1][-1]
                    g = windows[-1][0]
                    h = windows[-1][1]
                    i = windows[-1][-1]

                    dzdx[ax_y, ax_x] = np.nanmean([a,b,c,d,f,g,h,i])
                    
                    loc_string[ax_y, ax_x] = 'middle grid'
    return(dzdx,loc_string)
def nnmean_ax_ay(mat,ax_x,ax_y):
    '''mat =matrix has to be ar leest 3 x 3, x and y are turend'''
#     x = d3.x.values.copy()
#     y = d3.y.values.copy()
    ny,nx = mat.shape
   
    # gen = window3x3(bd)
    # windows_3x3 = np.asarray(list(gen))
    xmin_a,xmax_a,ymin_a,ymax_a = window3x3_a(mat)

    dzdx = np.empty((ny, nx)) #matrix for nanmean
  
    loc_string = np.empty((ny, nx), dtype="S25")

    windows = mat[xmin_a[ax_y]:xmax_a[ax_y],ymin_a[ax_x]:ymax_a[ax_x]]

    if ax_x == 0 and ax_y == 0:  # top left corner
        dzdx[ax_y, ax_x] = np.nanmean([windows[0][1] ,windows[1][0] , windows[1][1]])
        loc_string[ax_y, ax_x] = 'top left corner'
    elif ax_x == nx - 1 and ax_y == 0:  # top right corner
        dzdx[ax_y, ax_x] = np.nanmean([windows[1][1],windows[0][0],windows[1][0]])
               
        loc_string[ax_y, ax_x] = 'top right corner'

    elif ax_x == 0 and ax_y == ny - 1:  # bottom left corner
        dzdx[ax_y, ax_x] = np.nanmean([windows[1][1],windows[0][0],windows[0][1]])
                
        loc_string[ax_y, ax_x] = 'bottom left corner'

    elif ax_x == nx - 1 and ax_y == ny - 1:  # bottom right corner
        dzdx[ax_y, ax_x] = np.nanmean([windows[1][0],windows[0][0],windows[0][1]])
        loc_string[ax_y, ax_x] = 'bottom right corner'

       
    elif (ax_y == 0) and (ax_x != 0 and ax_x != nx - 1): # top boarder
        dzdx[ax_y, ax_x] = np.nanmean([windows[0][0],windows[0][-1],windows[1][0],windows[1][1],windows[1][2]])
        loc_string[ax_y, ax_x] = 'top boarder'

        
    elif ax_y == ny - 1 and (ax_x != 0 and ax_x != nx - 1):# bottom boarder
        dzdx[ax_y, ax_x] = np.nanmean([windows[1][0],windows[1][-1],windows[0][0],windows[0][1],windows[0][2]])
                    
        loc_string[ax_y, ax_x] = 'bottom boarder'

                # left boarder
    elif ax_x == 0 and (ax_y != 0 and ax_y != ny - 1):
        dzdx[ax_y, ax_x] = np.nanmean([windows[0][0],windows[1][-1],windows[2][0],windows[0][1],windows[2][1]])
        loc_string[ax_y, ax_x] = 'left boarder'

               
    elif ax_x == nx - 1 and (ax_y != 0 and ax_y != ny - 1): # right boarder
        dzdx[ax_y, ax_x] = np.nanmean([windows[0][0],windows[2][0],windows[0][1],windows[2][1],windows[1][0],])
                   
        loc_string[ax_y, ax_x] = 'right boarder'

            # middle grid
    else:
        a = windows[0][0]
        b = windows[0][1]
        c = windows[0][-1]
        d = windows[1][0]
        f = windows[1][-1]
        g = windows[-1][0]
        h = windows[-1][1]
        i = windows[-1][-1]

        dzdx[ax_y, ax_x] = np.nanmean([a,b,c,d,f,g,h,i])
                    
        loc_string[ax_y, ax_x] = 'middle grid'
        
    return(dzdx[ax_y, ax_x],loc_string)




def iter_circle_indices(tt):    
    x,y = tt.shape
    cy,cx=int((y-1)/2),int((x-1)/2)

    iter_x =[]
    iter_y =[]
    iter_x.append(cx)
    iter_y.append(cy)
    # goes r,u,l,d;
    # in x +1,
    # stay
    moves = np.linspace(1,x,x).astype('int')
    moves_r_u =moves[0::2]
    moves_l_d =moves[1::2]
    stop =0
    x_unique = np.arange(0,x)
    y_unique = np.arange(0,y)

    directions = ['r','u','d','l']
    for c,mv_max in enumerate(moves):
        d = c+1
        print(d)
        if d%2 ==1:
            if stop:
                break
            for mv in range(0,mv_max):

                #go right:
                if (iter_y[-1] +1) >(y-1):
                    print(' right upper boarder')
                    stop =1
                    break                
                else:
                    iter_x.append(iter_x[-1])
                    iter_y.append(iter_y[-1]+1)

            if stop:
                break
            for mv in range(0,mv_max):
                if ((iter_x[-1] -1) <0):
                    print('reached upper borader')
                    stop =1
                    break
                else:
                    iter_x.append(iter_x[-1] -1)
                    iter_y.append(iter_y[-1])

        else:
            if stop:
                break
            for mv in range(0,mv_max):
                #go left:
                if (iter_y[-1] -1) <0:
                    print('left borader')
                    stop =1
                    break
                else:
                    iter_x.append(iter_x[-1])
                    iter_y.append(iter_y[-1]-1)
            if stop:
                break

            for mv in range(0,mv_max):
                #go down:
                if (iter_x[-1] +1) >(x-1):
                    print('reached lower borader')
                    stop =1
                    break
                else:
                    iter_x.append(iter_x[-1]+1)
                    iter_y.append(iter_y[-1])

    iter_x_unique =np.unique(iter_x)
    iter_y_unique =np.unique(iter_y)
    dic_missing ={}
    if (len(iter_x_unique) != len(x_unique)):
        x_missing = np.setdiff1d(x_unique, iter_x_unique)
        print('x missing')
        dic_missing['x']=x_missing

    if (len(iter_y_unique) != len(y_unique)):
        y_missing = np.setdiff1d(y_unique, iter_y_unique)
        print('y missing')
        dic_missing['y']=y_missing

#     for key in dic_missing.keys():
#         if key == 'x':
#             for j,x_val in enumerate(dic_missing[key]):

#                 print (tt[x_val,])
#         else:
#             for j,y_val in enumerate(dic_missing[key]):

#                 print (tt[:,y_val])
    tupl_lst =tuple(zip(iter_x,iter_y))
    return(iter_x,iter_y,tupl_lst,dic_missing)
    
# for i in range 
