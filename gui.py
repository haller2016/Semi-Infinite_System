#ES3002 Project, by Eric Hall
#Last Updated, 3/6/2020

#Import important libraries
from appJar import gui #used for making GUI
import numpy as np #Used for making arrays 
from scipy import special #used for erf/erfc functions
import matplotlib as mpl #Used for plotting
import matplotlib.pyplot as plt #again used for plotting
import math #used for more advanced math functions
import time #used for keeping track of time

#Create GUI Object
app = gui('ES3002 Project, Eric Hall')

#Define important global variables (placeholders for now)
answers = 0 #Input responses
surf_c = 0 #Surface concentration (mol/m3)
initial_c = 0 #Initial Domain Concentration (mol/m3)
diff_co = 0 #Diffusion Coefficient (mol/m2*s)
domain_l = 0 #Characteristic Length of Domain (m)
z_range = 0 #Range of Z values created based on characterstic length
max_time = 0 #Max time calculated from semi-infinite rule
t_range = 0 #t-range calculated based on max_time variable
graph_num = 100 #used to restrict the size of the t_range and z_range variables
graph_speed = .01 #How quickly the graphs are produced
track_num = 0  #How many times the find stats tool has been used
point_plot = 0 #Plot of specific concentration at Z and T inputted
point_text = 0 #Annotation that gives stats on graph (flux, concentration)
line_plot = 0 #Plot of all z concentrations at a given time
flux_text = 0 #String version of Flux

#Function for checking validity of inputs (aka are they floats/ints instead of strings?)
def checkValid(box_name):
    value = app.getEntry(box_name) #Get value of box
    #If the value is blank, set it to the black star symbol (called wait in this library)
    if value == '': 
        app.setValidationEntry(box_name, state="wait")
    
    #If not...
    else:
        #Try getting the input value and converting it to a float. If it works, set to valid (green). If not, set to invalid (red). 
        try:
            good_conc = float(value) 
            app.setValidationEntry(box_name, 'valid')
        except ValueError:
            app.setValidationEntry(box_name, 'invalid')

#Running the plotting simulation
def runSim(button_name):
    #Pull in important globals used in this function
    global answers
    global surf_c
    global initial_c
    global diff_co
    global domain_l
    global z_range
    global max_time
    global t_range
    
    if button_name == 'Close':
        app.stop()
        plt.close()
    elif button_name == 'Run Sim':
        plt.clf()
        app.showSubWindow('Find Stats')
        try:
            answers = app.getAllEntries() #Get entries from input boxes
            surf_c = float(answers['Cs_Box']) #Surface Concentration
            initial_c = float(answers['C0_Box']) #Initial Domain Concentration
            diff_co = float(answers['Dab_Box']) #Diffusion Coeffiencent 
            domain_l = float(answers['length_Box']) #Characterstic Length of Domain
            #print(answers) 
            z_range = np.arange(0, domain_l, domain_l/graph_num) #Create range of Zs based on Characteristic Length and #of Plots wanted
            max_time = (domain_l**2) / diff_co #Find the maximum time based on the semi-infinite rule
            t_range = np.arange((max_time/1000), max_time, (max_time/graph_num)) #Create t-range based on number of graphs and max time

            #Calculate concentration for all times and zs!
            for t in t_range:
                phi = z_range/(2*math.sqrt(diff_co*t)) #Define phi (same for both uptake and downtake)
                if surf_c > initial_c:
                    theta = special.erfc(phi) #Caculate theta for all zs at the given time (uptake)
                    final_conc = (theta*(surf_c - initial_c)) + initial_c #Convert to concentration
                    plt.title('Semi-Infinite Uptake System') #Change Plot Title to Uptake
                else:
                    theta = special.erf(phi) #Calculate theta for all zs at the given time (release)
                    final_conc = (theta*(initial_c - surf_c)) + surf_c #Convert to concentration
                    plt.title('Semi-Infinite Release System') #Change plot title to Release

                #Plot Concentration vs Z at the given time, label axis
                plt.plot(z_range, final_conc, 'b', alpha = .35) 
                plt.xlabel('Z, [meters]')
                plt.ylabel('Concentration, [mol/m3]')

                plt.pause(.00000001)
                #print("Average Theta: " + str(sum(theta)/len(theta)))

        #Prints error in command line if an invalid input is given
        except ValueError:
            print('Invalid Input!')

#Find Concentration and Boundary Flux at a Certain Time    
def find_Conc():
    #Pull in global variables 
    global track_num 
    global point_plot
    global point_text
    global line_plot


    try:
        #Get specific entries from "Find Stats" Window
        stat_time = float(app.getEntry('Time')) #Specified Time
        stat_z = float(app.getEntry('Z-Coordinate')) #Specified Z-Coordinate
        line_phi = z_range/(2*math.sqrt(diff_co*stat_time)) #Create phi for all zs at this specific time (later shown in green)
        point_phi = stat_z/(2*math.sqrt(diff_co*stat_time)) #Create phi for the specific point (later shown as a black star)
        check_domain = math.sqrt(stat_time*diff_co) #Check to see if the time is valid within the domain (if not, line is changed to red!)

        #Uptake
        if surf_c > initial_c:
            line_theta = special.erfc(line_phi) #Calculate theta for all zs at the given time
            line_final_conc = (line_theta*(surf_c - initial_c)) + initial_c #Convert to concentration

            point_theta = special.erfc(point_phi) #Calculate theta for the specific point inputted
            point_final_conc = (point_theta*(surf_c - initial_c)) + initial_c #Convert to concentation
            
        #Release
        else:
            line_theta = special.erf(line_phi) #Calculate theta for all zs at the given time
            line_final_conc = (line_theta*(initial_c - surf_c)) + surf_c #Convert to concentation

            point_theta = special.erf(point_phi) #Calculate theta for the specific point inputted
            point_final_conc = (point_theta*(initial_c - surf_c)) + surf_c #Convert to concentation
            
        #If time is invalid, make graph red. Else, make green
        if check_domain > domain_l:
            graph_color = 'r'
        else:
            graph_color = 'g'
        
        dummy_num = .2*(surf_c-initial_c) #Number used to position text on graph

        current_flux = format((math.sqrt(diff_co/(math.pi*stat_time)) * (surf_c - initial_c)), '.3g') #Flux calculated for specfic time at Z = 0m

        string_conc = format(point_final_conc, '.3g') #Make concentration a string
        #If this is the first time getting stats...
        if track_num == 0:
            line_plot = plt.plot(z_range, line_final_conc, graph_color) #Plot line at specific t
            point_text = plt.annotate(f'Time: {stat_time}s\nZ-Coordinate: {stat_z}m\nConcentration: {string_conc}mol/m3\nMolar Boundary Flux, Na (z = 0): {current_flux}mol/m2*s', xy=(stat_z, point_final_conc), xytext=(domain_l*.4, surf_c - dummy_num), arrowprops=dict(facecolor='black')) #plot annotation with info
            point_plot = plt.plot(stat_z, point_final_conc, 'k*') #Plot specific concentration that was solved for earlier

        #If not the first time solving for stats...
        else:
            #Remove old annotations and graphs
            point_text.remove()
            plt.setp(point_plot, visible=False)
            plt.setp(line_plot, visible=False)
            
            #Plot New Graphs
            line_plot = plt.plot(z_range, line_final_conc, graph_color) #Plot line at specific t
            point_text = plt.annotate(f'Time: {stat_time}s\nZ-Coordinate: {stat_z}m\nConcentration: {string_conc}mol/m3\nMolar Boundary Flux, Na (z = 0): {current_flux}mol/m2*s', xy=(stat_z, point_final_conc), xytext=(domain_l*.4, surf_c - dummy_num), arrowprops=dict(facecolor='black')) #plot annotation with info
            point_plot = plt.plot(stat_z, point_final_conc, 'k*') #Plot specific concentration that was solved for earlier
        
        plt.pause(graph_speed) #Show graph without blocking
        track_num += 1 #Track the number of times stats has been solved for

    #Throw error is input is invald
    except ValueError:
        print('Invalid Input!')


#Make Labels/Boxes
app.addLabel('title', 'Semi-Infinite System Modeling', colspan=2) #Title of Main Window

app.addLabel('space_1', '', row=1) #Space for looks

#Add label, box, and set function for surface concentration
app.addLabel('Cs', 'Surface Concentration, [mol/m3]', row=2, column=0) 
app.addValidationEntry('Cs_Box', row=2, column=1)
app.setEntryChangeFunction('Cs_Box', checkValid)

#Add label, box, and set function for initial concentration
app.addLabel('C0', 'Initial Domain Concentration, [mol/m3]', row=3, column=0)
app.addValidationEntry('C0_Box', row=3, column=1)
app.setEntryChangeFunction('C0_Box', checkValid)

#Add label, box, and set function for diffusion coefficient
app.addLabel('Dab', 'Diffusion Coefficient, [m2/s]', row=4, column=0)
app.addValidationEntry('Dab_Box', row=4, column=1)
app.setEntryChangeFunction('Dab_Box', checkValid)

#Add label, box, and set function for domain length
app.addLabel('length', 'Domain Length, [m]', row=5, column=0)
app.addValidationEntry('length_Box', row=5, column=1)
app.setEntryChangeFunction('length_Box', checkValid)

app.addLabel('space_2', '', row=6) #another space for looks

#Add buttons and set their press function
app.addButtons(['Run Sim', 'Close'], runSim, row=7, colspan=2)

#Live-Stats Window title
app.startSubWindow('Find Stats', modal=False, blocking=False)

#Label for Find Stats Window
app.addLabel('live_stats', 'Fill out t and Z to find a Specific Concentration', row=0, colspan=2)

#Add label, box, and set function for time
app.addLabel('t', 'Time, [s]', row=1, column=0)
app.addValidationEntry('Time', row=1, column=1)
app.setEntryChangeFunction('Time', checkValid)

#Add label, box, and set function for Z-Coordinate
app.addLabel('z', 'Z-Coordinate, [m]', row=2, column=0)
app.addValidationEntry('Z-Coordinate', row=2, column=1)
app.setEntryChangeFunction('Z-Coordinate', checkValid)

#Add find stats button and set its function
app.addButton('Find Stats', find_Conc, row=3, colspan=2)

app.stopSubWindow()

#Start the GUI!
app.go()

