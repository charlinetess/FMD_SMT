####### This document plot the action map (prefered direction) generated from a data st generated by the attached file SMT.jl
# it plots the thing for 2 different trials 

# go to folder 
cd("Documents/FosterDayanMorris/StandardMemoryTest") # Define directory

# define function to compute place cells activity

# This function tells within wich index column is located x, used to take decision on which action to follow
function indice(Acum,x) # x number, Acum vector
    
    for i=1:length(Acum)
       if i==1
           if x<Acum[i] # if the random number generated is before the first 
                return i
            end
        else
            if Acum[i-1]<x<=Acum[i]
                return i
            end
        end
    end  
        
end
function  placecells(pos,centres,width)
# PLACECELLS(POSITION,CENTRES,WIDTH) calculates the activity of the place cells
#in the simulation. The returned vector F is of length N, where N is the number of place
#cells, and it contains the activity of each place cell given the simulated rat's current
#POSITION (a 2 element column vector). The activity of the place cells is modelled as a
#rate-of-fire (i.e. a scalar value) determined by a gaussian function. The CENTRES of the
#gaussian functions are an argument, and must be a 2 x N matrix containing each place
#cell's preferred location in 2D space. The WIDTH of the place cell fields must
#also be provided as a scalar value (all place cells are assumed to have the same
#width).
#
#The returned vector, F, must be a N element column vector.
    # calculate place cell activity
F = exp.(-sum((repeat(pos,1,size(centres,2))-centres).^2,dims=1)/(2*width^2))';
return F
end

# load packages 
using LinearAlgebra
using Statistics
using JLD2
using FileIO

# chose rat, day and trial you want to plot for.  
indexrat=8;
indextrial1=1;
indextrial2=18;
indexday=1;





# # initalize the value map variable
# vbegin = zeros(length(x),length(x));
# vend = zeros(length(x),length(x));
widthplacecells=[0.4]*100;
discountfactor=[0.98];

for l=1:length(widthplacecells)
for k=1:length(discountfactor)
    let widthplacecells=[0.4]*100,discountfactor=[0.98], ubegin, vbegin, uend, vend, steps, vbegin,x,y,x2,y2,vend,Wbegin,Wend,R,r, parameters,data, centres, width;


# load data 
rats=load("experiment_$(widthplacecells[l])_$(discountfactor[k]).jld2");
parameters=rats["parameters"];
featuresexperiment=rats["features"]; # conteins the fields numberofdays numberoftrials and numberofrats
data=rats["data"];
#data[indexrat][indexday][indextrial] contains fields :  trajectory, latency, searchpreference,actionmap,valuemap,TDerror,PCcentres 

r=parameters[:r]; # platform radius
R=parameters[:R];
angles=parameters[:angles];
temperature=parameters[:temperature]; # parameter in the exponential 



# establish the grid of points in the pool
steps=20;
x=[-R+(steps)*(k-1) for k=1:(2*R/steps+1)];
y=zeros(1,length(x));
transpose!(y,x);
x2=x;
y2=y;


# initalize the vector map variables
ubegin = zeros(length(x),length(x));
vbegin = zeros(length(x),length(x));
uend = zeros(length(x),length(x));
vend = zeros(length(x),length(x));


zbegin=data[indexrat][indexday].day[indextrial1].actionmap;
zbegin=0*data[indexrat][indexday].day[indextrial1].actionmap;

zend=data[indexrat][indexday].day[indextrial2].actionmap;


centres=parameters[:centres];
width=widthplacecells[l];
width=parameters[:σPC];

# for each place point in the grid, calculate the vector of preferred action direction
for i = 1:length(x)
    for j = 1:length(x)
        # make sure the point is in the pool
        if sqrt((x[i]^2+x[j]^2)) < R
            # determine the place cell activity at this point
            F = placecells([x[i],y[j]],centres,width);     
            #  Compute action cell activity    
            actactioncellend=transpose(zend)*F; 
            #  Compute action cell activity    
            actactioncellbegin=transpose(zbegin)*F;
            if maximum(actactioncellbegin)>=100 
                actactioncellbegin=100*actactioncellbegin./maximum(actactioncellbegin); 
                elseif maximum(actactioncellend)>=100
                actactioncellend=100*actactioncellend./maximum(actactioncellend); 
            end
             
            # Compute probability distribution : 
            Pactioncellbegin=exp.(temperature.*actactioncellbegin)./sum(exp.(temperature.*actactioncellbegin)); 
            Pactioncellend=exp.(temperature.*actactioncellend)./sum(exp.(temperature.*actactioncellend)); 

            # Compute summed probability distribution:
           SumPactioncellbegin=[sum(Pactioncellbegin[1:k]) for k=1:length(Pactioncellbegin)    ]
           SumPactioncellend=[sum(Pactioncellend[1:k]) for k=1:length(Pactioncellend)    ]
           # Generate uniform number between 0 and 1 :
           xbegin=rand();
           xend=rand();
           # now chose action: 
           indexactionbegin=indice(SumPactioncellbegin,xbegin); # Chose which action between the 8 possibilities
           indexactionend=indice(SumPactioncellend,xend); # Chose which action between the 8 possibilities
           argdecisionbegin=angles[indexactionbegin]; # compute the corresponding angle 
           argdecisionend=angles[indexactionend]; # compute the corresponding angle 
           newdirbegin=[cos(argdecisionbegin) sin(argdecisionbegin)];
           newdirend=[cos(argdecisionend) sin(argdecisionend)];
           # store the result in u and v
            ubegin[i,j] = 10*newdirbegin[1];
            vbegin[i,j] = 10*newdirbegin[2];
            uend[i,j] = 10*newdirend[1];
            vend[i,j] = 10*newdirend[2];
        else
            #x[i] = NaN;
            #y[j] = NaN;
            ubegin[i,j]= NaN;
            vbegin[i,j] = NaN;
            uend[i,j]= NaN;
            vend[i,j] = NaN;
        end
    end
end

theta=0:pi/50:(2*pi+pi/50); # to plot circles 
# Plot value function : 
 using PyPlot


clf()
ioff()


# # create the figure 
# fig = figure("Test plot action map rat $(indexrat)",figsize=(10,5));
# suptitle("action map rat $(indexrat), after 1 trial $(indextrial1), trial $(indextrial2)")

# subplot(121)
# plot(R*cos.(theta),R*sin.(theta),"k-")
# plot(data[indexrat][indexday].platformposition[1].+r*cos.(theta),data[indexrat][indexday].platformposition[2].+r*sin.(theta),"m-")
# quiver(x,y,ubegin,vbegin,color="b");
# xlabel("X Position (cm)");
# ylabel("Y Position (cm)");

# ax=gca() 
# ax[:set_axis_off]()


# subplot(122)
# plot(R*cos.(theta),R*sin.(theta),"k-")
# plot(data[indexrat][indexday].platformposition[1].+r*cos.(theta),data[indexrat][indexday].platformposition[2].+r*sin.(theta),"m-")
# quiver(x,y,uend,vend,color="b");
# xlabel("X Position (cm)");
# ylabel("Y Position (cm)");
# ax=gca() 
# ax[:set_axis_off]()
# show()




fig = figure("Action Map",figsize=(6,6));
        # plot circle 
        plot(R*cos.(theta),R*sin.(theta),color="darkgrey",lw=1)
plot(data[indexrat][indexday].platformposition[1].+r*cos.(theta),data[indexrat][indexday].platformposition[2].+r*sin.(theta),color="slateblue")
quiver(x,y,uend,vend,color="darkred");
# quiver(x,y,ubegin,vbegin,color="darkred");

xlabel("X Position (cm)");
ylabel("Y Position (cm)");
ax=gca() 
ax[:set_axis_off]()
        #gca()[:grid](false);
        #gca()[:view_init](20.0,0.0)



# lets try this 
# cbar=fig[:colorbar]#(ax=axe,image)
#cbar.ax.tick_params(labelsize=20) 


show()


        savefig("Action_$(parameters[:σPC])_$(parameters[:γ]).png")

      end  # end scope variables 
    end 
  end 


#savefig("Actionmap$(rats.parameters).png")

