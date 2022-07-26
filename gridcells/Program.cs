﻿// See https://aka.ms/new-console-template for more information
using gridcells;

var spatial = new SpatialNavigation();


var conv = spatial.Conv(45);


Console.WriteLine(conv.Item1);
Console.WriteLine(conv.Item2);

spatial.RandomNavigation(3);
spatial.Plot();


var grid = new Grid();
var simulation = new Simulation(grid, spatial.txx, spatial.tyy);
simulation.run();

Console.WriteLine("End");
