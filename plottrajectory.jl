# load data 
using JLD2
using FileIO
# rats=load("/Users/pmxct2/Documents/FosterDayanMorris/Sublime/experimentDMP.jld2");

rats=load("experiment_40.0_0.98.jld2");
parameters=rats["parameters"];
featuresexperiment=rats["features"];# contains the fields numberofdays numberoftrials and numberofrats
data=rats["data"];

#data[indexrat][indexday][indextrial] contains fields :  trajectory, latency, searchpreference,actionmap,valuemap,TDerror,PCcentres 


r=5; # platform radius
R=100



# chose rat 
indexrat=8;
# chose Day
indexday=1;
# chose trial
indextrial=1;



argument=0:pi/50:2pi+pi/50;
xplat=r*cos.(argument);
yplat=r*sin.(argument);
xmaze=R*cos.(argument);
ymaze=R*sin.(argument);


#Declare a figure object 
using PyPlot
ioff()
fig = figure("trajectory",figsize=(9,9))





ax1 = gca() 
# grid("on") # Create a grid on the axis
#title("trial $(indextrial1)") # Give the most recent axis a title
ax1[:set_ylim]([-101,101])
ax1[:set_xlim]([-101,101])
xlabel("X")
ylabel("Y")


# Plot place cells 
#scatter(parameters[:centres][1,:],parameters[:centres][2,:],s=3)
# plot platform
plot(data[indexrat][indexday].platformposition[1].+xplat,data[indexrat][indexday].platformposition[2] .+ yplat,color="red")

# Plot circle
plot(xmaze,ymaze,color="mediumseagreen",lw=2)


# Plot trajectory 
plot(data[indexrat][indexday].day[indextrial].trajectory[:,1],data[indexrat][indexday].day[indextrial].trajectory[:,2],color="darkred", lw=2)
ax1.set_axis_off()

show()


