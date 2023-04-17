
function benchmark_systems_ff(ground)
    [
        ("cyclic 7",  Groebner.cyclicn(7, ground=ground)), 
        ("cyclic 8",  Groebner.cyclicn(8, ground=ground)), 
        ("cyclic 9",  Groebner.cyclicn(9, ground=ground)), 
        ("katsura 10",Groebner.katsuran(10, ground=ground)), 
        ("katsura 11",Groebner.katsuran(11, ground=ground)), 
        ("katsura 12",Groebner.katsuran(12, ground=ground)), 
        ("eco 11",    Groebner.eco11(ground=ground)),         
        ("eco 12",    Groebner.eco12(ground=ground)),         
        ("eco 13",    Groebner.eco13(ground=ground)),         
        ("noon 7",    Groebner.noonn(7, ground=ground)),  
        ("noon 8",    Groebner.noonn(8, ground=ground)),  
        ("noon 9",    Groebner.noonn(9, ground=ground)),   
        ("henrion 5", Groebner.henrion5(ground=ground)),  
        ("henrion 6", Groebner.henrion6(ground=ground)),  
        ("henrion 7", Groebner.henrion7(ground=ground)),  
        ("reimer 6",  Groebner.reimern(6, ground=ground)), 
        ("reimer 7",  Groebner.reimern(7, ground=ground)), 
        ("reimer 8",  Groebner.reimern(8, ground=ground)), 
    ] 
end

function benchmark_systems_qq(ground)
    [
        ("cyclic 7",  Groebner.cyclicn(7, ground=ground)), 
        ("cyclic 8",  Groebner.cyclicn(8, ground=ground)), 
        ("katsura 9",Groebner.katsuran(9, ground=ground)), 
        ("katsura 10",Groebner.katsuran(10, ground=ground)), 
        ("eco 10",    Groebner.eco10(ground=ground)),         
        ("eco 11",    Groebner.eco11(ground=ground)),         
        ("noon 8",    Groebner.noonn(8, ground=ground)),  
        ("noon 9",    Groebner.noonn(9, ground=ground)),   
        ("henrion 6", Groebner.henrion6(ground=ground)),  
    ] 
end
