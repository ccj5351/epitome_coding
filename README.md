# Epitome Transform Coding: Towards Joint Compression of a Set of Images

## Overall Framework of Our Approach

<img src="files/workflow-1.png" style="width: 25%;"/>
<img src="files/workflow-2.png" style="width: 25%;"/>

- Learning the epitome $E$ of a collection of images ${I_i }$, and doing the reconstruction via $E$ and the associated transform map ${\Psi_i}$.

- Then the bitstream of entropy-encoded epitome, transform maps, and residuals, can be transmitted with bandwidth saving and economic storage.

- Last, the transmitted bitstream will be decoded for the final rendering.

## Epitome Learning and Reconstruction

### What is the Epitome?
- An epitome is a condensed image of parameters that specifies a generative model of patches taken from
input images.
- It contains high-order statistics of the texture and shape properties of the input images.
- A epitome in size $W_e \times H_e$ can be viewed as a mixture of $W_e \cdot H_e$ Gaussians.
- The epitome is learned so that if small patches are sampled in an unordered fashion from it, they will have nearly the same appearance as patches sampled from the original input image.

### Epitome Encoding

### Transform Map Encoding

### Residual Processing and Encoding
