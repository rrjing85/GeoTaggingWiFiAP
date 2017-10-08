import math

from math import radians, cos, sin, asin, sqrt

from math import atan, atan2, sin, cos, tan, pi, sqrt, radians, degrees

def linear_median(arr):
	arr = sorted(arr)
	if len(arr) % 2 != 0:
		return arr[(len(arr))/2]
	else:
		return float(arr[len(arr)/2] + arr[(len(arr ) - 1)/2])/2

def candMedian(dataPoints):
	#Calculate the first candidate median as the geometric mean 
	tempX = 0.0
	tempY = 0.0
	tempZ = 0.0
	
	for i in range(0,len(dataPoints)): 
		tempX += dataPoints[i][0] 
		tempY += dataPoints[i][1]
		tempZ += dataPoints[i][2]
	return [tempX/len(dataPoints), tempY/len(dataPoints), tempZ/len(dataPoints)]

def numersum(testMedian,dataPoint):
	# Provides the denominator of the weiszfeld algorithm depending on whether you are adjusting the candidate x or y 
	return 1/math.sqrt((testMedian[0]-dataPoint[0])**2 + (testMedian[1]-dataPoint[1])**2 + (testMedian[2]-dataPoint[2])**2)

def denomsum(testMedian, dataPoints):
	# Provides the denominator of the weiszfeld algorithm 
	temp = 0.0
	for i in range(0,len(dataPoints)):
		temp += 1/math.sqrt((testMedian[0] - dataPoints[i][0])**2 + (testMedian[1] - dataPoints[i][1])**2 + (testMedian[2] - dataPoints[i][2])**2) 
	return temp

def objfunc(testMedian, dataPoints):
	# This function calculates the sum of linear distances from the current candidate median to all points # in the data set, as such it is the objective function we are minimising.
	temp = 0.0
	for i in range(0,len(dataPoints)):
		temp += math.sqrt((testMedian[0]-dataPoints[i][0])**2 + (testMedian[1]-dataPoints[i][1])**2 + (testMedian[2]-dataPoints[i][2])**2) 
	return temp

def objfunc_distances(distances):
	return sum(distances)

def denomsum_distances(distances):
	temp = 0.0
	for distance in distances:
		temp += 1/distance
	return temp

def getdistances(testMedian, dataPoints):
	distances = []
	for i in range(0,len(dataPoints)):
		distances.append(math.sqrt((testMedian[0]-dataPoints[i][0])**2 + (testMedian[1]-dataPoints[i][1])**2 + (testMedian[2]-dataPoints[i][2])**2)) 
	return distances

#Calculates the median of an array of 3D points	
def getmedian(dataPoints):
	testMedian = candMedian(dataPoints)
	numIter = 50
	last_objfunc = 0.0
	for x in range(0,numIter):
		#print objfunc(testMedian,dataPoints)
		distances = getdistances(testMedian, dataPoints)
		current_objfunc = objfunc_distances(distances)
		if math.fabs(last_objfunc - current_objfunc) < 1:
			break
		last_objfunc =  current_objfunc
		denom = denomsum_distances(distances) 
		nextx = 0.0
		nexty = 0.0
		nextz = 0.0
		
		for y in range(0,len(dataPoints)):
			numer = 1/distances[y]#numersum(testMedian,dataPoints[y])
			nextx += (dataPoints[y][0] * numer)/denom 
			nexty += (dataPoints[y][1] * numer)/denom
			nextz += (dataPoints[y][2] * numer)/denom
		testMedian = [nextx,nexty,nextz]
	
	return testMedian

#Calculate median absolute deviation

def get_MAD_for_points_with_median(points, median):
	diffs = []
	for element in points:
		diffs.append(math.sqrt((element[0] - median[0])**2 + (element[1] - median[1])**2 + (element[2] - median[2])**2))
	mad_python = linear_median(diffs)
	return mad_python

def getMAD(points):
	if len(points) < 2:
		return -1
	median = getmedian(points)
	return get_MAD_for_points_with_median(points, median)

def getMAD_median(arr):
	mad_python = -1
	median_value = arr[0]
	if len(arr) > 1:
		median_value = getmedian(arr)
		diffs = []
		for element in arr:
			diffs.append(math.sqrt((element[0] - median_value[0])**2 + (element[1] - median_value[1])**2 + (element[2] - median_value[2])**2))
		mad_python = linear_median(diffs)
	return (mad_python, median_value)


def getMean(arr):
	sum_x = 0
	sum_y = 0
	sum_z = 0
	for element in arr:
		sum_x = sum_x + element[0]
		sum_y = sum_y + element[1]
		sum_z = sum_z + element[2]
	div = len(arr)
	return (float(sum_x)/div, float(sum_y)/div, float(sum_z)/div)

def getSTD(arr):
	mean = getMean(arr)
	diffs = []
	for element in arr:
		diffs.append((element[0] - mean[0])**2 + (element[1] - mean[1])**2 + (element[2] - mean[2])**2)
	return math.sqrt(float(sum(diffs))/len(diffs))



def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r


f = 1/298.257223
R = 6378137
e = 8.1819190842622e-2

def llaToECEF(coords):
    
	coords = [radians(val) for val in coords]
	(lat, lon, h) = coords
	lambds = atan((1-f)**2*tan(lat))
	rs = sqrt(R**2/(1+(1/(1-f)**2-1)*sin(lambds)**2))
	
	x = rs*cos(lambds)*cos(lon) + h*cos(lat)*cos(lon)
	y = rs*cos(lambds)*sin(lon) + h*cos(lat)*cos(lon)
	z = rs*sin(lambds) + h*sin(lat)
	
	return (x,y,z)

def ECEFTolla(coords):
	(x, y, z) = coords
	b = sqrt(pow(R,2) * (1-pow(e,2)))
	ep = sqrt((pow(R,2)-pow(b,2))/pow(b,2))
	p = sqrt(pow(x,2)+pow(y,2))
	th = atan2(R*z, b*p)
	lon = atan2(y, x)
	lat = atan2((z+ep*ep*b*pow(sin(th),3)), (p-e*e*R*pow(cos(th),3)))
	n = R/sqrt(1-e*e*pow(sin(lat),2))
	alt = p/cos(lat)-n
	lat = (lat*180)/pi
	lon = (lon*180)/pi
	
	return (lat, lon, alt)
	
# ecef = llaToECEF(radians(12.5774271), radians(55.7847283), radians(0))
# print ecef
# print ECEFTolla(ecef[0], ecef[1], ecef[2])


from pylab import plt
import numpy as np
import networkx as nx

def partition_to_draw(partition_obj):
    if type(partition_obj == 'dict'):
        partition = np.array([partition_obj[i] for i in range(len(partition_obj))])
    
    count = 0.
    list_nodes = []
    color      = []
    for com in set(partition_obj):
        count = count + 1.

        this_com_nodes = []
        for i in range(len(partition_obj)):      
            if partition_obj[i] == com:
                this_com_nodes.append(i)

        list_nodes.extend(this_com_nodes)

        color.extend([count] * len(this_com_nodes))

    color = np.array(color)
    
    return list_nodes, color
    


def draw_partitioned_graph(G, partition_obj, layout=None, labels=None,layout_type='spring', 
               node_size=70, node_alpha=0.7, cmap=plt.get_cmap('jet'),
               node_text_size=12,
               edge_color='blue', edge_alpha=0.5, edge_tickness=1,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    # if a premade layout haven't been passed, create a new one
    if not layout:
        if graph_type == 'spring':
            layout=nx.spring_layout(G)
        elif graph_type == 'spectral':
            layout=nx.spectral_layout(G)
        elif graph_type == 'random':
            layout=nx.random_layout(G)
        else:
            layout=nx.shell_layout(G)

    # prepare the partition list noeds and colors

    list_nodes, node_color = partition_to_draw(partition_obj)
      
    # draw graph
    nx.draw_networkx_nodes(G,layout,list_nodes,node_size=node_size, 
                           alpha=node_alpha, node_color=node_color, cmap = cmap)
    nx.draw_networkx_edges(G,layout,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    #nx.draw_networkx_labels(G, layout,font_size=node_text_size,
    #                        font_family=text_font)

    if labels is None:
        labels = range(len(G))

    edge_labels = dict(zip(G, labels))
    #nx.draw_networkx_edge_labels(G, layout, edge_labels=edge_labels, 
    #                            label_pos=edge_text_pos)

    # show graph

    plt.axis('off')
    plt.xlim(0,1)
    plt.ylim(0,1)








