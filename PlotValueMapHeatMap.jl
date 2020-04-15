#### define functions 


# # define function to compute place cells activity
# function  placecells(pos,cent,width)
# #
# # PLACECELLS(POSITION,CENTRES,WIDTH) calculates the activity of the place cells
# #in the simulation. The returned vector F is of length N, where N is the number of place
# #cells, and it contains the activity of each place cell given the simulated rat's current
# #POSITION (a 2 element column vector). The activity of the place cells is modelled as a
# #rate-of-fire (i.e. a scalar value) determined by a gaussian function. The CENTRES of the
# #gaussian functions are an argument, and must be a 2 x N matrix containing each place
# #cell's preferred location in 2D space. The WIDTH of the place cell fields must
# #also be provided as a scalar value (all place cells are assumed to have the same
# #width).
# #
# #The returned vector, F, must be a N element column vector.
#     # calculate place cell activity

# F = exp.(-sum((repeat(pos,1,size(cent,2))-cent).^2,dims=1)./(2*width.^2));
# Fbis=zeros(length(F),1)
# transpose!(Fbis,F)
# return Fbis
# end


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

# load data 

# go to folder 
cd("Documents/FosterDayanMorris/StandardMemoryTest") # Define directory

# load every package 
using LinearAlgebra
using Statistics
using JLD2
using FileIO

using PyPlot

# chose rat, trials and day 

indexrat=1;

indextrial1=1;
indextrial2=20;

indexday=1;

# intialise arrays 

# establish the grid of points in the pool
global steps
steps=2;

# # initalize the value map variable
# vbegin = zeros(length(x),length(x));
# vend = zeros(length(x),length(x));

widthplacecells=[0.1 0.3 0.4 0.5 0.98];
discountfactor=[0.1 0.4 0.6 0.9 0.7];

for l=1:length(widthplacecells)
for k=1:length(discountfactor)
    let widthplacecells=[0.1 0.3 0.4 0.5 0.98]*100, width=widthplacecells[l], discountfactor=[0.1 0.4 0.6 0.9 0.7], vbegin,x,y,x2,y2,vend,Wbegin,Wend,R,r;

        # load data 
        rats=load("experiment_$(widthplacecells[l])_$(discountfactor[k]).jld2");

        parameters=rats["parameters"];
        featuresexperiment=rats["features"];

        data=rats["data"];

        #data[indexrat][indexday][indextrial] contains fields :  trajectory, latency, searchpreference,actionmap,valuemap,TDerror,PCcentres 


        r=parameters[:r]; # platform radius
        R=parameters[:R]; # maze radius 
        theta=parameters[:angles];

        x=[-R+(steps)*(k-1) for k=1:(2*R/steps+1)];
        y=zeros(1,length(x));
        transpose!(y,x);
        x2=x;
        y2=y;
        # initalize the value map variable
        vbegin = zeros(length(x),length(x));
        vend = zeros(length(x),length(x));

        Wbegin=data[indexrat][indexday].day[indextrial1].valuemap;
        Wend=data[indexrat][indexday].day[indextrial2].valuemap;

        # if PCplast==0
        width=widthplacecells[l];
        # widthsend=0.30*100*ones(1,size(featuresexperiment[4],2));
        centres=parameters[:centres];
        # centresend=featuresexperiment[4];
        # elseif PCplast==1
        # widthsbegin=data[indexrat][indexday].day[indextrial1].PCwidths;
        # widthsend=data[indexrat][indexday].day[indextrial2].PCwidths;

        # centresbegin=data[indexrat][indexday].day[indextrial1].PCcentres;
        # centresend=data[indexrat][indexday].day[indextrial2].PCcentres;

        # end

        # for each place point in the grid, calculate the critic value
        # for i = 1:length(x)

        #     for j = 1:length(x)

        #         # make sure the point is in the pool
        #         if sqrt((x[i]^2+y[j]^2)) < R
                
        #             # determine the place cell activity at this point
        #             F = placecells([x[i],y[j]],centres,width)     
        #             # Fend = placecells([x[i],y[j]],centresend,widthsend)       
        #             # determine the actor activity
        #             vbegin[j,i] = dot(Wbegin,F)[1];
        #             vend[j,i] = dot(Wend,F)[1];
        #         else
        #             vbegin[j,i] = NaN;
        #             vend[j,i] = NaN;
        #         end
        #     end
        # end
        # this is commented as it does not show the correct value map for the initialisation, that is 0 everywhere, but already after one trial update. However culd be interesting to compare this for different values of gamma 






        for i = 1:length(x)

            for j = 1:length(x)

                # make sure the point is in the pool
                if sqrt((x[i]^2+y[j]^2)) < R
                
                    # determine the place cell activity at this point
                    F= placecells([x[i],y[j]],centres,width)       
                    # determine the actor activity
                    vbegin[j,i] = 0;
                    vend[j,i] = dot(Wend,F)[1];
                else
                    vbegin[j,i] = NaN;
                    vend[j,i] = NaN;
                end
            end
        end


        # plot value function :


        # # Define color :
        # using Colors
        # using PyPlot,Colors

        # homemadecolor=ColorMap("A", [RGB(144/255,238/255,144/255),RGB(60/255,179/255,113/255),RGB(0,0.5,0)])
        # homemadecolorbis=ColorMap("B", [RGB(32/255,178/255,170/255),RGB(60/255,179/255,113/255),RGB(0,0.5,0)])
        # homemadecoloragain=ColorMap("C", [RGB(102/255,205/255,170/255),RGB(60/255,179/255,113/255),RGB(0,0.5,0)])
        # homemadecoloragain=ColorMap("C", [RGB(224/255,255/255,255/255),RGB(60/255,179/255,113/255),RGB(0,0.5,0)])
        # in the end using this colormap it is hard to see, the defaut one is better - at least for this 

        clf()
        ioff()
        # create the figure 
        fig = figure("Value function",figsize=(6,12));
        #suptitle("value map")


        ax=subplot(211)
        # show the value map
        #s1 = surf(x,y,vbegin,cmap=ColorMap("jet"));

        pcolormesh(x,x,vbegin)#,cmap=homemadecoloragain)#,aspect_ratio=1)


        # plot circle 
        plot(R*cos.(theta),R*sin.(theta),ls="--",color=[169/255,169/255,169/255],zorder=1)
        # plot platform
        plot(data[indexrat][indexday].platformposition[1].+r*cos.(theta),data[indexrat][indexday].platformposition[2].+r*sin.(theta),color=[250/255,128/255,114/255],zorder=2)



        xlabel("X Position (cm)");
        ylabel("Y Position (cm)");
        title("Before Learning")
        ax[:set_axis_off]()


        ax=subplot(212)


        # plot circle 
        plot(R*cos.(theta),R*sin.(theta),ls="--",color=[169/255,169/255,169/255],zorder=1)
        # plot platform
        plot(data[indexrat][indexday].platformposition[1].+r*cos.(theta),data[indexrat][indexday].platformposition[2].+r*sin.(theta),color=[250/255,128/255,114/255],zorder=2)


        pcolormesh(x2,x2,vend)#,cmap=homemadecoloragain)#,aspect_ratio=1)

        #fig[:canvas][:draw]()

        xlabel("X Position (cm)");
        ylabel("Y Position (cm)");
        title("After Learning")
        colorbar()
        #ax=gca() 
        ax[:set_axis_off]()
        #gca()[:grid](false);
        #gca()[:view_init](20.0,0.0)


# lets try this 
cbar=fig[:colorbar]#(ax=axe,image)
#cbar.ax.tick_params(labelsize=20) 

        savefig("Values_$(parameters[:σPC])_$(parameters[:γ]).png")

end # end scope variables 
end 
end 

show()



