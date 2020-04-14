# This code implements a code proposed by Foster Morris and Dayan 2000, Hippocampus. It implements the first model present in their article, in which a reinforcement learning agent is trained to learn a position of the platform. 

# it creates a .jld2 file containing the parmaters, the feature experiments and the data. Data contains trajectory, latency, actor weights and critic weights for every trial.

# This has been coded in Julia and works in julia version 1.1. 
# IT first calls the required packages. If you do not have them installed, open Pkg by ctrl + ] and type add NameOfPackage
# It then defines some functions that will be used within the loop 
# The section " parameters" contains all the parameters that can be played on to study the effect on the learning and on the dynamics
# Feature experiments defines the experiment features that is the number fo independant runs, the number of trials per day and the number of days 
# first define the directory you want to work in, and then just afyer having defined the parameters change the name of the file in whic you want to save the data on 



using LinearAlgebra
using Statistics
using JLD2
using FileIO

#cd("Documents/FosterDayanMorris/StandardMemoryTest/")

###################################################################################
###################################################################################
########################                                ###########################
########################       DEFINE FUNCTIONS         ###########################
########################                                ###########################
########################                                ##########################
########################                                ##########################
###################################################################################
###################################################################################


# generate place cell actibity uniformely distributed throughout the space
# this is actually unecessary - random distribution is enough 
#The algorithm places n points, of which the kth point is put at distance sqrt(k-1/2) from the boundary (index begins with k=1), and with polar angle 2*pi*k/phi^2 where phi is the golden ratio. Exception: the last alpha*sqrt(n) points are placed on the outer boundary of the circle, and the polar radius of other points is scaled to account for that. This computation of the polar radius is done in the function radius.
#function  radius(k,n,b) # k index of point on the boundary, n number of total points, b number of boundary points
#    if k>n-b
#        r = 1;            # put on the boundary
#    else
#        r = sqrt(k-1/2)/sqrt(n-(b+1)/2);     # computation of radius of the different points 
#    end
#end


# sunflower seed arrangement :
#function sunflower(n, R, alpha)   # n number of centers,
#    # alpha is indicating how much one cares about the evenness of boundary , chose 2 to have nice trade off
#    # R is the radius of the circle in cm
#
#    b = round(alpha*sqrt(n));      # number of boundary points
#    phi = (sqrt(5)+1)/2;           # golden ratio
#    
#    r=[R*radius(k,n,b) for k=1:n];
#    theta=[2*pi*k/phi^2 for k=1:n];
#    
#    # scatter(r.*cos.(theta), r.*sin.(theta));#, marker='o', "m");
#    X=r.*cos.(theta); 
#    Y=r.*sin.(theta);
#    return hcat(X, Y)
#end


# Define activity as a function of position 

# Compute the activity of the place cells
function  placecells(pos,centres,width)
#
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



# Calculate reward as a function of position 
function reward(x,y,xp,yp,parameters) # x,y position of the rat and xp,yp position of the platform, r radius of the platform
    if (x-xp)^2+(y-yp)^2<= parameters[:r]^2 # if the rat is in the platform
        R=parameters[:rewardgoal];
    else # else 
        R=0;
    end 
end


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




#########################################################################
#############          LOOP       ###############  ######################
#########################################################################

function STM(parameters,featuresexperiment,NameOfFile)  
  experiment=[];
    for indexrat=1:featuresexperiment[:numberofrats]
            ##########  ##########  ##########  ##########   ########## 
        ##########  ##########  START EXPERIMENT  ##########  ##########  
            ##########  ##########  ##########  ##########   ########## 
            let currentexperiment, actorweights, criticweights # define scope of variables 
               criticweights=zeros(parameters[:NPC],1);
               actorweights=zeros(parameters[:NPC],parameters[:NA]);    
               currentexperiment=[];
                   for indexday=1:featuresexperiment[:numberofdays] # start experiment 
                       let indexplatform, xp, yp, currentday # start their scope 
                           # Everyday the location of the platform changes
                           # Chose platform :
                           indexplatform=rand(1:8); # take ith platform 
                           xp=parameters[:Xplatform][indexplatform];
                           yp=parameters[:Yplatform][indexplatform]; 
                           #xp=40;
                           #yp=40;
                                 ##########  ##########  ##########  ##########  
                           ##########  ##########  START DAY ##########  ##########  
                               ##########  ##########  ##########  ##########  
                           currentday=[];
                           for indextrial=1:featuresexperiment[:numberoftrials] ##########  
                               
                               ## Chose starting position :
                               #        
                               #indexstart=rand(1:4); # take indexstart-th starting position : chose    randomnly between 4 possibilities 1 East 2 North 3 West 4 South
                               #position=[Xstart[indexstart] Ystart[indexstart]];
                               
                               # start scope of variables: 
                               let indexstart, positionstart, currentposition, re,k,t,historyX,historyY,    TDerrors,arg,timeout,prevdir,dir
                                   # Chose starting position :     
                                   indexstart=rand(1:4); # take indexstart-th starting position : chose     randomnly between 4 possibilities 1 East 2 North 3 West 4 South
                                   positionstart=[parameters[:Xstart][indexstart] parameters[:Ystart][indexstart]];
                                   currentposition=positionstart;
                                   
                                   # Initialize reward 
                                   re=0;
                                   
                                   # Initialise index to save the trajectory and the values 
                                   k=1;
                                   # initialise time 
                                   t=parameters[:times][k];
                                   historyX=Float64[];
                                   historyY=Float64[];
                                   #valuemap=Float64[];
                                   TDerrors=Float64[];
                                   arg=0;        
                                   timeout=0;        
                                   prevdir=[0 0];    
                                   ##########  ##########  ##########  ##########   ########## 
                                   ##########  ##########  START TRIAL ##########  ##########  
                                   ##########  ##########  ##########  ##########   ########## 
                                   
                                       while t<=parameters[:T] && re==0
                           
                                               if t==parameters[:T]
                                                   X=xp;
                                                   Y=yp;
                                                   currentposition=[X Y];
                                                   timeout=1; # if we have to put the rat on the platform     then we dont reinforce the actor but only the critic
                                               end
                                           # Store former position to be able to draw trajectory
                                           push!(historyX,currentposition[1]); 
                                           push!(historyY,currentposition[2]);
                                                # compute new activity of pace cells :
                                           # actplacecell=place_activity(   position[1],position[2],Xplacecell,Yplacecell,σ); # this    function is wrong 
                                           actplacecell=placecells([currentposition[1],currentposition[2]],   parameters[:centres],parameters[:σPC]);
                                       
                                           ### Compute Critic ###
                                           C=dot(criticweights,actplacecell); # current estimation of the     future discounted reward 
                                           
                                           ####### Take decision and move to new position : ########
                                           # Compute the activity of action cells 
                           
                                           #  Compute action cell activity    
                                           actactioncell=transpose(actorweights)*actplacecell; # careful z    contains place cells in rows and action cells in column 
                                               if maximum(actactioncell)>=100
                                                   actactioncell=100*actactioncell./maximum(actactioncell); 
                                               end
                                           # Compute probability distribution : 
                                           Pactioncell=exp.(parameters[:temperature]*actactioncell)./sum(exp.(parameters[:temperature]*actactioncell)); 
                                           # Compute summed probability distribution:
                                           #SumPactioncell=cumul(Pactioncell);
                                           SumPactioncell=[sum(Pactioncell[1:k]) for k=1:length(Pactioncell)    ]
                       
                                           # Compute summed probability distribution:
                                           # SumPactioncell=cumul(Pactioncell); # other possibility 
                                           # Generate uniform number between 0 and 1 :
                                           x=rand();
                       
                                           # now chose action: 
                                           indexaction=indice(SumPactioncell,x); # Chose which action     between the 8 possibilities
                                           argdecision=parameters[:angles][indexaction]; # compute the coreesponding     angle 
                                           newdir=[cos(argdecision) sin(argdecision)];
                                           dir=(newdir./(1.0+parameters[:momentum]).+parameters[:momentum].*prevdir./(1.0+parameters[:momentum]));     # smooth trajectory to avoid sharp angles
                                               if !(norm(dir)==0)
                                                   #global dir
                                                   dir=dir./norm(dir); # normalize so we control the exact    speed of the rat
                                               end
                                           formerposition=currentposition;
                                           # Compute new position : 
                                           currentposition=currentposition.+parameters[:dt].*parameters[:speed].*dir; 
                                           if currentposition[1]^2+currentposition[2]^2>=parameters[:R]^2 # if we are outside of circle, move a lil bit 
                                               currentposition = (currentposition./norm(currentposition))*(parameters[:R]     - parameters[:R]/50);
                                               if !(norm(currentposition-formerposition)==0)
                                                   dir=(currentposition-formerposition)./norm(    currentposition-formerposition);
                                               else
                                                   dir=[0 0]
                                               end
                                           end
                                           prevdir=dir;
                                           ###  Compute reward ### 
                                           re=reward(currentposition[1],currentposition[2],xp,yp,parameters);                            
                                           # compute new activity of pace cells :
                                           # actplacecell=place_activity(   position[1],position[2],Xplacecell,Yplacecell,σ);
                                           actplacecell=placecells([currentposition[1],currentposition[2]],   parameters[:centres],parameters[:σPC]);
                       
                                               if re==1 # if we are on the platform 
                                                  ###  Compute error ###
                                                   Cnext=0;
                                               else 
                                                   Cnext=dot(criticweights,actplacecell);# new estimation     of the future discounted reward 
                                               end 
                                           
                                       
                                           #### Compute error  ####
                       
                                           err=re+parameters[:γ]*Cnext-C[1];
                       
                                           push!(TDerrors,err);
                                       
                                       
                                           ######### Compute new weights : ########
                                           
                                           # Actor weights :
                                               if timeout==0
                                                   G=zeros(8,1); # creating a matrix to select the row to     update
                                                   G[indexaction]=1; # indicating which row is to update 
                                                   # weights between action cells and place cells only    reinforced when the rats actually found the platform
                                                   #    z[:,indexaction]=z[:,indexaction]+Z.*err.*actplacecell;   # only the weights between place cells and the action    taken are updated
                                                   actorweights=actorweights+parameters[:actorLR].*err.*actplacecell*    transpose(G);       
                                               end
                                           
                                           # Critic weights : 
                                           criticweights=criticweights+parameters[:criticLR].*err.*actplacecell;
                                           k=k+1; # counting steps
                                           t=parameters[:times][k]; # counting time
                                       
                                       ##################################################            
                                       end
                                       ########## ##########  END TRIAL ########## ##########             
                                   push!(historyX,currentposition[1]); # Store the last position visited 
                                   push!(historyY,currentposition[2]);
                                   # push!(valuemap,w)
                                               
                                              ############### SAVING THE THINGS IN THE DIFFERENT CLASS    ################
                                   ## in creating a new trial type one should write Trial(Trajectory,     latency, actionmap) # action map atm is just z, then    it will be improved adding a new attribute being value map 
                                   
                                   currenttrial=(trajectory=hcat(historyX,historyY),latency=t,actionmap=actorweights,valuemap=    criticweights,TDerror=TDerrors); # Creating the current trial with all    its fields
                                   push!(currentday,currenttrial);# Storing it in the current day 
                               end # end scope of variables 
                           ##################################################     
                           end 
                           ########## ##########  END DAY ########## ##########
                            ########## ##########  END DAY ########## ##########
                           dayc=(day=currentday, platformposition=[xp, yp]);
                           push!(currentexperiment,dayc);       
                       end # End scope of variable       
                   ##################################################     
                   end 
                   ########## ##########  END EXPERIMENT ########## ##########
               push!(experiment,currentexperiment);
           end # end scope 
    ##################################################     
    end # end all the rats 
save(NameOfFile, "parameters",parameters,"features",featuresexperiment,"data",experiment);
end # end function






widthplacecells=[0.1 0.3 0.4 0.5 0.98];
discountfactor=[0.1 0.4 0.6 0.9 0.7];

for i=1:4
for j=1:4
	let parameters,featuresexperiment,NameOfFile, widthplacecells=[0.05 0.1 0.2 0.7], discountfactor=[0.1 0.4 0.6 0.9]
###################################################################################
################################### PARAMETERS  ###################################
###################################################################################

# Creating the circle and the place cells:
R= 100; # Radius of the circle in cm
r=5;# Radius of the platform  in cm

# Motion characteristic 
dt=0.1; # timestep in s 
speed=30; # speed of the rat in cm.s-1
# Different possible directions 
angles=[-3*pi/4, -2*pi/4, -pi/4, 0, pi/4, 2*pi/4, 3*pi/4, pi];


# Trial characteristic :
T=120; # maximal duration of a trial in seconds

# Place cells 
NPC=493; # number of place cells 
#centres=sunflower(N,R,2);
#Xplacecell=centres[:,1]; # absciss place cells  
#Yplacecell=centres[:,2]; # y place cells 


# Place cell : method used by Blake richards 
# initialize the centres of the place cells by random unifrom sampling across the pool
arguments= rand(1,NPC)*2*pi; # generate random angles 
radii= sqrt.(rand(1,NPC))*R; # generate random radius 
centres= [cos.(arguments).*radii; sin.(arguments).*radii];  # gather 


σPC=widthplacecells[j]*100;

#σPC=0.30*100; # variability of place cell activity, in centimeters
ampρPC=1;

# Action cells : 
NA=8; # number of action cells 


# Potential positions of the platform : 
Xplatform=[0.3,0,-0.3,0,0.5,-0.5,0.5,-0.5].*R; # in cm
Yplatform=[0,0.3,0,-0.3,0.5,0.5,-0.5,-0.5].*R;# in cm

# Potential Starting positions of the rat :
Xstart=[0.95,0,-0.95,0].*R; # East, North, West, South
Ystart=[0,0.95,0,-0.95].*R;

# Define number of rats, number of days and numbers of trials per day
#global numberofdays, numberoftrials, numberofrats
numberofdays=2;
numberofrats=30;
numberoftrials=20;

#global times
times=collect(0:dt:T+dt); # time vector 


# Parameter that regulate the choice between former angle and new angle 
#global momentum
momentum=1.0;

temperature=2; # in reality inverse temperature, if high more exploitation, low more exploration 
# Learning variables : 
#global γ, actorLR, criticLR, criticLR
γ=discountfactor[i];
#γ=0.98; # Discount factor.  they dont precise the value  
actorLR=0.1; # actor learning rate
criticLR=0.01; # critic learning rate

rewardgoal=1; # value of reward given when finding platform 

#########################################################################
#############        Create input     ###############  ######################
#########################################################################


parameters=Dict(:momentum=>momentum,:γ=>γ,:actorLR=>actorLR,:criticLR=>criticLR,:centres=>centres,:R=>R,:r=>r,:speed=>speed,:angles=>angles,:NPC=>NPC,:NA=>NA,:σPC=>σPC,:ampρPC=>ampρPC,:Xstart=>Xstart,:Ystart=>Ystart,:dt=>dt,:T=>T,:times=>times,:Xplatform=>Xplatform,:Yplatform=>Yplatform,:rewardgoal=>rewardgoal,:temperature=>temperature);
featuresexperiment=Dict(:numberofrats=>numberofrats, :numberofdays=>numberofdays, :numberoftrials=>numberoftrials);





NameOfFile="experiment_$(σPC)_$(γ).jld2";



@time  STM(parameters,featuresexperiment,NameOfFile)

end # end let 
end # end loop over i 

end # end loop over j 
