# load data 

# go to folder 
# cd("Documents/FosterDayanMorris/StandardMemoryTest") # Define directory

# load every package 
using LinearAlgebra
using Statistics
using JLD2
using FileIO

using PyPlot



# for i=1:length(widthplacecells)
# 	for j=1:length(discountfactor)
# 		let widthplacecells=[0.4]*100,discountfactor=[0.1]

# load data 
rats=load("experiment_40.0_0.98.jld2");
parameters=rats["parameters"];
featuresexperiment=rats["features"];# contains the fields numberofdays numberoftrials and numberofrats
data=rats["data"];
#data[indexrat][indexday][indextrial] contains fields :  trajectory, latency, searchpreference,actionmap,valuemap,TDerror,PCcentres 

# load parameters so we dont have to call the structure all the time 
r=parameters[:r]; # platform radius
R=parameters[:R]; # maze radius 


# Define number of rats, number of days and numbers of trials per day
numberofdays=featuresexperiment[:numberofdays];
numberofdays=1
numberofrats=featuresexperiment[:numberofrats];
numberoftrials=featuresexperiment[:numberoftrials]; 

# Computing the mean of all latencies 
latencies=[mean([data[n][div(k+numberoftrials-1,numberoftrials)].day[rem(numberoftrials-1+k,numberoftrials)+1].latency for n in 1:numberofrats]) for k in 1:numberoftrials*numberofdays ];




clf()
ioff()
fig = figure("Test plot latencies",figsize=(9,9))
#ax = fig[:add_subplot](1,1,1)

   
global labels_code=[]

for k=0:(numberofdays-1)

	# Calculate standard deviation 
	#err=[std([rats.experiment[n].day[k].trial[i].Latency for n in 1:numberofrats]; corrected=false) for i in 1:numberoftrials] ;

	# Calculate the lower value for the error bar : 
	uppererror = [std([data[n][k+1].day[i].latency for n in 1:numberofrats]; corrected=false) for i in 1:numberoftrials]./sqrt(numberofrats) ;
	lowererror = [std([data[n][k+1].day[i].latency for n in 1:numberofrats]; corrected=false) for i in 1:numberoftrials]./sqrt(numberofrats) ;

	errs=[lowererror,uppererror];

	PyPlot.plot(k*numberoftrials.+(0:numberoftrials-1), [mean([data[n][k+1].day[i].latency for n in 1:numberofrats]) for i in 1:numberoftrials ], marker="None",linestyle="-",color="darkgreen",label="Base Plot")


	PyPlot.errorbar(k*numberoftrials.+(0:numberoftrials-1),[mean([data[n][k+1].day[i].latency for n in 1:numberofrats]) for i in 1:numberoftrials ],yerr=errs,fmt="o",color="k")
	global labels_code=vcat(labels_code,collect(1:numberoftrials))
	println(labels_code)

end 
ax=gca() 

mx = matplotlib.ticker.MultipleLocator(1) # Define interval of minor ticks
ax.xaxis.set_major_locator(mx) # Set interval of minor ticks

ax.spines["top"].set_color("none")
ax.spines["right"].set_color("none")
ax.spines["bottom"].set_visible("False")
ax.spines["left"].set_visible("False")

xmin, xmax = ax.get_xlim() 
ymin, ymax = ax.get_ylim()
# get width and height of axes object to compute 
# matching arrowhead length and width
dps = fig.dpi_scale_trans.inverted()
bbox = ax.get_window_extent().transformed(dps)
width, height = bbox.width, bbox.height
# manual arrowhead width and length
hw = 1/20*(ymax-ymin) 
hl = 1/20*(xmax-xmin)
lw = 1 # axis line width
ohg = 0.3 # arrow overhang
# compute matching arrowhead length and width
yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width 
yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height

ax.arrow(xmin, ymin, xmax-xmin, 0.,length_includes_head= "True", fc="k", ec="k", lw = lw,head_width=hw, head_length=hl, overhang = ohg,  clip_on = "False") 

ax.arrow(xmin, ymin, 0., ymax-ymin,length_includes_head= "True", fc="k", ec="k", lw = lw, head_width=yhw, head_length=yhl, overhang = ohg,  clip_on = "False")

xlabel("Trials ")#, fontsize=18);
ylabel("Latencies")#, fontsize=18)

#labels = [item.get_text() for item in ax.get_xticklabels()]
#labels[1] = labels_code
labels_code=vcat("",labels_code)
ax.set_xticklabels(labels_code)

SMALL_SIZE = 10
MEDIUM_SIZE = 20
BIGGER_SIZE = 30

plt.rc("font", size=SMALL_SIZE)          # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)    # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title



savefig("Latencies.png")

# 		end # end let 
# 	end 

# end 

#arrowed_spines(fig, ax)
show()


# 	
# 	                              ,,          ,,            ,,
# 	                            `7MM   mm     db          `7MM
# 	                              MM   MM                   MM
# 	`7MMpMMMb.pMMMb.`7MM  `7MM    MM mmMMmm `7MM `7MMpdMAo. MM  .gP"Ya
# 	  MM    MM    MM  MM    MM    MM   MM     MM   MM   `Wb MM ,M'   Yb
# 	  MM    MM    MM  MM    MM    MM   MM     MM   MM    M8 MM 8M""""""
# 	  MM    MM    MM  MM    MM    MM   MM     MM   MM   ,AP MM YM.    ,
# 	.JMML  JMML  JMML.`Mbod"YML..JMML. `Mbmo.JMML. MMbmmd'.JMML.`Mbmmd'
# 	                                               MM
# 	                                             .JMML.



# discountfactors=[0.1 0.3 0.4 0.5 0.6 0.9 0.98];
widths=[5.0 10.0 30.0 40.0 50.0 70.0];
# 		let widthplacecells=[0.4]*100,discountfactor=[0.1]

lat=[]; # create empty latency array 
errors=[];

# for discountfactor in discountfactors 
for width in widths 
	# load data 
	# rats=load("experiment_40.0_$(discountfactor).jld2");
	rats=load("experiment_$(width)_0.98.jld2");
	parameters=rats["parameters"];
	featuresexperiment=rats["features"];# contains the fields numberofdays numberoftrials and numberofrats
	data=rats["data"];
	#data[indexrat][indexday][indextrial] contains fields :  trajectory, latency, searchpreference,actionmap,valuemap,TDerror,PCcentres 

	# Define number of rats, number of days and numbers of trials per day
	numberofdays=featuresexperiment[:numberofdays];
	numberofdays=1
	numberofrats=featuresexperiment[:numberofrats];
	numberoftrials=featuresexperiment[:numberoftrials]; 
	# Computing the mean of all latencies 
	latencies=[mean([data[n][div(k+numberoftrials-1,numberoftrials)].day[rem(numberoftrials-1+k,numberoftrials)+1].latency for n in 1:numberofrats]) for k in 1:numberoftrials*numberofdays ];

	uppererror =[std([data[n][div(k+numberoftrials-1,numberoftrials)].day[rem(numberoftrials-1+k,numberoftrials)+1].latency for n in 1:numberofrats])./sqrt(numberofrats) for k in 1:numberoftrials*numberofdays ]

	#uppererror = [[std([data[n][k+1].day[i].latency for n in 1:numberofrats]; corrected=false) for i in 1:numberoftrials]./sqrt(numberofrats) for k=1:numberofdays];

	lat=push!(lat,latencies) # lat[indexdiscount] to access data 
	errors=push!(errors,uppererror); # errors[indexdiscount][indexdays] is the way to access the data 
end 





colors=["firebrick","peru", "darkorange","forestgreen","darkturquoise","royalblue","darkviolet"];
# or darkturquoise before royalblue or dhotpink at the end 



clf()
ioff()
fig = figure("Test plot latencies",figsize=(9,9))
#ax = fig[:add_subplot](1,1,1)

   
global plotlegend=[] # need to create am array of plots for the legend 

for k=1:length(lat)
	p1,=plot(collect(1:numberoftrials), lat[k], marker="None",linestyle="-",color=colors[k],label=latexstring("\$\\gamma=$(discountfactors[k])\$")) # if we forget the comma afte the assignment for some reason it doesnt work when plotting the legend # because the plot assignment is a list, we must assign the first element of the lsit to p1 
	errs=vcat(transpose(errors[k]),transpose(errors[k]))
	PyPlot.errorbar(collect(1:numberoftrials),lat[k],yerr=errs,fmt="o",color=colors[k])
	push!(plotlegend,p1)
end 



ax=gca() 

mx = matplotlib.ticker.MultipleLocator(1) # Define interval of minor ticks
ax.xaxis.set_major_locator(mx) # Set interval of minor ticks

ax.spines["top"].set_color("none")
ax.spines["right"].set_color("none")
ax.spines["bottom"].set_visible("False")
ax.spines["left"].set_visible("False")

ax[:set_ylim]([0,125])
ax[:set_xlim]([0,21])
xmin, xmax = ax.get_xlim() 
ymin, ymax = ax.get_ylim()
# get width and height of axes object to compute 
# matching arrowhead length and width
dps = fig.dpi_scale_trans.inverted()
bbox = ax.get_window_extent().transformed(dps)
width, height = bbox.width, bbox.height
# manual arrowhead width and length
hw = 1/20*(ymax-ymin) 
hl = 1/20*(xmax-xmin)
lw = 1 # axis line width
ohg = 0.3 # arrow overhang
# compute matching arrowhead length and width
yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width 
yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height

ax.arrow(xmin, ymin, xmax-xmin, 0.,length_includes_head= "True", fc="k", ec="k", lw = lw,head_width=hw, head_length=hl, overhang = ohg,  clip_on = "False") 

ax.arrow(xmin, ymin, 0., ymax-ymin,length_includes_head= "True", fc="k", ec="k", lw = lw, head_width=yhw, head_length=yhl, overhang = ohg,  clip_on = "False")

xlabel("Trials ")#, fontsize=18);
ylabel("Latencies")#, fontsize=18)

#labels = [item.get_text() for item in ax.get_xticklabels()]
#labels[1] = labels_code
labels_code=vcat(collect(1:numberoftrials))
labels_code=vcat("","",labels_code)
ax.set_xticklabels(labels_code)

SMALL_SIZE = 10
MEDIUM_SIZE = 20
BIGGER_SIZE = 30

plt.rc("font", size=SMALL_SIZE)          # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)    # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

# categories = [latexstring("\$\\gamma=$(discountfactors[i])\$") for i=1:length(discountfactors)]
categories = [latexstring("\$\\sigma=$(widths[i])\$") for i=1:length(widths)]

# leg3 = legend(vcat(p5,p[1:length(TimePower)],p5,p[length(TimePower)+1,end]),vcat([latexstring("\$\\sigma=$(sigma)\$")] , categories , [latexstring("\$\\sigma=$(sigma2)\$")] , categories),loc=2, ncol=2) # Two columns, vertical group labels

# leg3 = legend(vcat(p5,pright,p5,pleft),vcat([latexstring("\$\\sigma=$(sigma)\$")] , categories , [latexstring("\$\\sigma=$(sigma2)\$")] , categories),loc=2, ncol=2) # Two columns, vertical group labels
p2, = plot([0], marker="None",
           linestyle="None", label="dummy-tophead")
# leg4 = legend(vcat(p2,plotlegend),vcat([latexstring("\$\\gamma\$")] , categories),loc=1)# ncol=1) # Two columns, vertical group labels 
leg4 = legend(vcat(p2,plotlegend),vcat([latexstring("\$\\sigma\$")] , categories),loc=1)
 legend(plotlegend,categories,loc=1)


show()




