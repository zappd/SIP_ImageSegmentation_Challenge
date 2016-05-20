# SIP Image Segmentation Challenge

This repo contains a complete implementation of the paper [Spectral Embedding and Min-Cut for Image
Segmentation](http://www.cs.utoronto.ca/~strider/publications/minCutSpec.pdf) for **Hyperspectral Images**. It includes an implementation written in C to perform spectral embedding for hyperspectral images, as well as perform the necessary graph cuts.

You will need to link with [FEAST](http://www.ecs.umass.edu/~polizzi/feast/download.htm), [LAPACKE](http://www.netlib.org/lapack/#_software), and [vlfeat](https://github.com/vlfeat/vlfeat) to run.

For more information, see Implementation_Notes.pdf. Also note that for larger images, finding eigenvalues can take quite a bit of time.