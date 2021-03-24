using StatsPlots
using Plots.PlotMeasures

function valTQ(a)
    y = zeros(tmax-tmin+1)
    for i in a
        s, f, v = i
        for t in s-tmin+1:f-tmin
            y[t] = v
        end
    end
    return y
end

function showTQ(a; fc=:blue, xr=0, wh=[600, 250], fs=8, ti="TQ", bm=10px)
    ticklabel = string.(collect(tmin:tmax))
    y = valTQ(a)
    bar(y,
#        xlabel = "time", 
#        ylabel = "value",
        title = ti,
        fillcolor = fc,
        grid = false,
        lw = 0,
        bar_position = :stack,
        bar_width=1,
        legend = false,
        xrotation = xr,
        tickfontsize = fs,
        size = wh,
        bottom_margin=bm,
        xticks=(1:tmax-tmin+1, ticklabel)
    )    
end

function showPairTQ(a, b; xr=0, wh=[600, 250], fs=8, ti="TQ pair", 
    fc = [:royalblue1 :tomato1 :purple4], bm=10px)
    ticklabel = string.(collect(tmin:tmax))
    y1 = valTQ(a)
    y2 = valTQ(b)
    m = min.(y1,y2)
    z1 = y1-m
    z2 = y2-m
    groupedbar([z1 z2 m], 
        title = ti,
        bar_position = :stack, 
        fillcolor = fc, 
        label = ["a" "b" "ab"],
        grid = false,
        lw = 0,
        legend = false,
        xrotation = xr,
        tickfontsize = fs,
        size = wh,
        bottom_margin=bm,
        bar_width=1,
        xticks=(1:tmax-tmin+1, ticklabel)
    )
end
