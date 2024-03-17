
point_coordinate = {'K':[0.33333333333,0.3333333333,0],'M':[0,0.5,0],'G':[0,0,0]}

def band_points (points):
    point_list =[]
    for point in points:
        coordinate = symmetry_coordinate(point)
        point_list.append(coordinate)
    return point_list

def symmetry_coordinate(point):
    point = point.upper()
    coordinate = point_coordinate[point]
    return coordinate

# def band_input(points,num_points):
#     coordinates = band_points(points)
#     parameter = []
#     for j,i in enumerate(coordinates):
#         input_line = {'x':str(i[0]),'y':str(i[1]),'z':str(i[2]),'number':str(num_points),'label':f" ! {points[j].upper()}"}
#         parameter.append(input_line)
#     return parameter


def band_input(path,points,num_points):
    parameter = []
    for k,i in enumerate(path):
        for j in points:
            if i.upper() == j[0].upper():
                if k == len(path)-1:
                    input_line = {'x':str(j[1]),'y':str(j[2]),'z':str(j[3]),'number':str(1),'label':f" ! {j[0].upper()}"}
                else:
                    input_line = {'x':str(j[1]),'y':str(j[2]),'z':str(j[3]),'number':str(num_points),'label':f" ! {j[0].upper()}"}
                parameter.append(input_line)
    return parameter
