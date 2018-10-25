#!/usr/bin/env python3

import sys
import os
import subprocess


filter_script_mlx = """<!DOCTYPE FilterScript>
<FilterScript>
 <filter name="Compute normals for point sets">
  <Param description="Neighbour num" type="RichInt" name="K" value="10"/>
  <Param description="Smooth Iteration" type="RichInt" name="smoothIter" value="0"/>
  <Param description="Flip normals w.r.t. viewpoint" type="RichBool" name="flipFlag" value="true"/>
  <Param z="0" description="Viewpoint Pos." y="0" type="RichPoint3f" name="viewPos" x="0"/>
 </filter>
 <filter name="Normalize Vertex Normals"/>
 <filter name="Marching Cubes (RIMLS)">
  <Param name="FilterScale" type="RichFloat" description="MLS - Filter scale" value="2"/>
  <Param name="ProjectionAccuracy" type="RichFloat" description="Projection - Accuracy (adv)" value="0.0001"/>
  <Param name="MaxProjectionIters" type="RichInt" description="Projection - Max iterations (adv)" value="15"/>
  <Param name="SigmaN" type="RichFloat" description="MLS - Sharpness" value="0.75"/>
  <Param name="MaxRefittingIters" type="RichInt" description="MLS - Max fitting iterations" value="3"/>
  <Param name="Resolution" type="RichInt" description="Grid Resolution" value="200"/>
 </filter>
 <filter name="Change the current layer">
  <Param description="Mesh" type="RichMesh" name="mesh" value="0"/>
 </filter>
 <filter name="Delete Current Mesh"/>
 <filter name="Merge Close Vertices">
  <Param max="2.7852" description="Merging distance" min="0" type="RichAbsPerc" name="Threshold" value="0.027852"/>
 </filter>
 <filter name="Remove Duplicate Faces"/>
 <filter name="Remove Duplicate Vertices"/>
 <filter name="Remove Unreferenced Vertices"/>
 <filter name="Discrete Curvatures">
  <Param enum_cardinality="4" value="0" enum_val2="RMS Curvature" enum_val3="ABS Curvature" type="RichEnum" enum_val0="Mean Curvature" enum_val1="Gaussian Curvature" description="Type:" name="CurvatureType"/>
 </filter>
</FilterScript>"""


# <filter name="Marching Cubes (APSS)">
#  <Param description="MLS - Filter scale" type="RichFloat" name="FilterScale" value="2"/>
#  <Param description="Projection - Accuracy (adv)" type="RichFloat" name="ProjectionAccuracy" value="0.0001"/>
#  <Param description="Projection - Max iterations (adv)" type="RichInt" name="MaxProjectionIters" value="15"/>
#  <Param description="MLS - Spherical parameter" type="RichFloat" name="SphericalParameter" value="1"/>
#  <Param description="Accurate normals" type="RichBool" name="AccurateNormal" value="true"/>
#  <Param description="Grid Resolution" type="RichInt" name="Resolution" value="200"/>
# </filter>





cwd = os.getcwd()

def create_tmp_filter_file(filename='filter_file_tmp.mlx'):
    curwordir = cwd.split('/')
    index = curwordir.index("MATLAB")
    tempdir = '/'.join(curwordir[0:index+1])
    try:
        os.mkdir(tempdir + "/tmp")
    except OSError as e:
        print(sys.stderr, "Exception creating folder for meshes: " + str(e))

    with open(tempdir + '/tmp/' + filename, 'w') as f:
        f.write(filter_script_mlx)
    return tempdir + '/tmp/' + filename



def run_meshlab(in_file, out_file,
                 filter_script_path=create_tmp_filter_file()):
    # # Add input mesh
    command = "meshlab.meshlabserver -i " + in_file
    # Add the output filename and output flags
    command += " -o " + out_file + " -m vn vq"
    # Add the output options
    # Add the filter script
    command += " -s " + filter_script_path
    # Execute command
    print("Going to execute: " + command)
    output = subprocess.call(command, shell=True)
    print()
    print("Done:")
    print(in_file, " > ", out_file, ": ", str(output))






if __name__ == '__main__':
    # if len(sys.argv) < 3:
    #     print("Usage:")
    #     print(sys.argv[0] + " /path/to/input_mesh num_iterations")
    #     print("For example, reduce 10 times:")
    #     print(sys.argv[0] + " /home/myuser/mymesh.dae 10")
    #     exit(0)

    in_mesh_path = sys.argv[1]
    in_mesh = in_mesh_path.split('/')[-1]
    filedir = '/'.join(in_mesh_path.split('/')[0:-1])
    out_mesh = sys.argv[2]
    out_mesh_path = filedir + "/" + out_mesh

    print("Directory: ", filedir)
    print("Input mesh: ", in_mesh_path, " (filename: " , in_mesh, ")")
    print("Output mesh: ", out_mesh_path, " (filename: " , out_mesh, ")")

    # try:
    #     os.mkdir(tmp_folder_name)
    # except OSError as e:
    #     print( sys.stderr, "Exception creating folder for meshes: ", str(e))


    run_meshlab(in_mesh_path, out_mesh_path)


    # in_mesh = sys.argv[1]
    # print("in_mesh ", in_mesh)
    # filename = in_mesh.split('/')[-1]
    # print("filename: ", filename)
    #
    #
    # folder_name = filename.replace('.', '_')
    # print("folder_name: ", folder_name)
    # tmp_folder_name = cwd + '/tmp/' + folder_name + '_meshes/'
    # print("tmp_folder_name: ", tmp_folder_name)
    #
    # print("Input mesh: " + in_mesh + " (filename: " + filename + ")")
    # print("Output folder: " + tmp_folder_name)
    #
    # try:
    #     os.mkdir(tmp_folder_name)
    # except OSError as e:
    #     print(sys.stderr, "Exception creating folder for meshes: " + str(e))
