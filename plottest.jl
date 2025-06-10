using CSV
using DataFrames
using Plots; gr()

species = ["bS","bS2","bS3","bS4","bS5","bS6","bS7","bS8"]

obs = [2.285e-6,2.144e-9,9.714e-8,1.142e-8]

colors = [:black :cyan :magenta :green]
obshape = [:circle]


for i = 2:length(species)+1
  file = "ab/$(species[i-1]).ab"
  f = open(file)
  lines = readlines(f)
  timesteps = length(lines)
  temp = Array{Float64,timesteps}
  if i == 2 
    global plotarray = zeros(timesteps,1+length(species))
#    global obsarray = zeros(timesteps,1+length(species))
  end
  for j = 1:length(lines)
    if j != 1 
      # On the first run, populate the time column
      if i == 2 
        plotarray[j,i-1] = log10(parse(Float64,split(lines[j])[1]))
#        obsarray[j,i-1]  = log10(parse(Float64,split(lines[j])[1])) 
      end
      plotarray[j,i]   = log10(parse(Float64,split(lines[j])[2]))
#      obsarray[j,i]   = log10(obs[i-1])
    end
  end
end

plot(
     plotarray[2:end,1],
     plotarray[2:end,2:end],
     label=species, 
#     color = colors,
     xlabel = "log10(time (yr))",
     ylabel = "log10(n(x)/nH)",
     linewidth = 4,
     legend=:none
    )

#=
plot!(
      obsarray[2:end,1],
      obsarray[2:end,2:end],
      color = colors,
      line = :dot,
      linewidth = 4,
      label=species,
      legend=:best
     )
=#
