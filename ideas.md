# Headroom analysis and next ideas

This note asks where an anchor compressor could still win on images from `gen_image.py`. It uses three kinds of evidence. A calculator gives the ideal cost of any coding order on the generator's background (`lab/gauss_orders.py`). The reference coder's cost per pixel was measured on synthetic images whose objects are known (`lab/gen_labeled.py`, `lab/costmap.c`, `lab/where_bits.py`). A literature search covered anchor selection, coding order and structural side information. The probes behind every measured number are listed in `lab/notes.md`, section 8. The ideas in the second half are untested, and their gains are estimates.

In short, coding order is worth nothing on the background, even with finite stencils, so no choice of anchor pixels can help there. The reference coder is about 3 KB per image above the ideal. A third of that is in the background, where the shared predictor is pulled off course by the objects, and the rest sits at outlines, lines, rings and striped discs. Edge detection helps only if it is used to describe objects compactly or to adapt the model to local edges, and both routes have to pay for their code twice.

## The ideal cost of the background does not depend on the order

Outside objects, `gen_image.py` makes a stationary Gaussian field: three power-law fields with fixed exponents and weights, multiplied by a scale of 2.15 to 3.0 grey levels per unit, plus white sensor noise and rounding. For such a field the best predictor from any set of known pixels is linear, and its error follows from the spectrum. A pixel then costs 0.5 log2(2πe·mse) bits, which matches the exact entropy of a rounded Gaussian to 0.001 bpp in the range that matters. On a synthetic background, one least-squares fit with a radius-7 stencil reaches 1.127 bpp, against 1.122 predicted.

Ideal cost at a scale of 2.5, close to `image.bin`:

| order | stencil radius | bpp |
|---|---|---|
| raster | 2 / 3 / 4 / 7 / 15 px | 1.2279 / 1.2107 / 1.2067 / 1.2035 / 1.2024 |
| raster, whole causal past (Szegő limit) | | 1.2021 |
| quadtree from a step-2 / 4 / 8 lattice | 6 px | 1.2027 / 1.2029 / 1.2035 |
| quadtree from a step-2 lattice | 3 / 2 px | 1.2052 / 1.2104 |
| lattice of step 16 / 4, then raster reading ahead onto the anchors | 7 px | 1.2048 / 1.2027 |
| every other pixel known | 15 px | 1.1105 |

Every order with a stencil of radius 6 or more lands within 0.003 bpp of the limit, less than 100 bytes per image. The levels of a quadtree cost very different amounts, from 1.56 bpp for anchors on a step-4 lattice to 1.12 bpp for the last edge midpoints, and they average out to the raster figure. The chain rule holds here even for finite stencils. Hierarchies win only when the stencil is tiny. At radius 2 the quadtree is 0.018 bpp better than raster. Memon, Neuhoff and Shende found the same ranking for small contexts in 2000.

## Where the reference coder loses bits

Two synthetic images, seeds 2 and 3 with archives of 39,517 and 36,156 bytes, were coded with the reference model and the cost of every pixel was grouped by what the pixel shows. The ideal for background is the Szegő limit at that image's scale. Inside opaque shapes the field is smoother, and its own limit is about 0.78 bpp. Outlines, lines, rings and stripes are fixed by a few parameters each, so their ideal is taken as the background's.

| pixels | share | bytes per image | bytes above ideal |
|---|---|---|---|
| background | 63.9% | 24,455 | 961 |
| inside opaque shapes | 14.8% | 3,849 | 257 |
| outlines of discs, ellipses and rectangles | 5.0% | 2,578 | 731 |
| lines | 3.8% | 1,795 | 424 |
| rings | 3.0% | 1,285 | 210 |
| striped discs | 9.4% | 3,873 | 400 |
| all | | 37,837 | 2,983 |

Outline pixels with a step of at least one grey level cost about 1.75 bpp, against about 1.2 for background. Spread over the 3-pixel band, that is about 1.6 extra bits per pixel of outline length.

The background loss has two causes. On a background-only image the coder's predictor has 3.4 to 4.6% more squared error than the ideal, which is the price of learning online and costs about 0.022 bpp. In real images its error on clean background rises to 12.6% above the ideal, because the same weights also learn object interiors, which want very different weights. With background weights the inside of a shape costs 0.055 bpp more than it should, and the reverse costs 0.10 bpp. This drift costs about 550 bytes per image. Probes ruled out three other suspects: weighting large errors down, solving the least-squares system more exactly, and correcting the Gaussian's variance for rounding all changed the archive by less than 65 bytes. The Gaussian input on its own is far from ideal on clean background, 1.244 bpp against 1.122, and the counters and the mixer recover most of that, to 1.170. Separately, the outer 8 pixels of the image cost 0.07 to 0.1 bpp more than the rest of the background, roughly 100 to 200 bytes per image.

## Why anchors placed by edges cannot pay

Positions are expensive. The best coders of sparse masks need about 4.4 to 8 bits per point at densities of 1 to 10% (Mohideen, Peter and Weickert 2021). The earlier probe in `lab/notes.md`, section 4, handed the decoder 4,000 anchors placed by edge strength for free, and the other pixels got cheaper by 385 bits in total, less than 0.1 bit per anchor.

An anchor pixel also carries the wrong kind of information. On the background no set of pixels helps, as shown above. At an outline the extra cost is spread along the whole length. An anchor fixes the edge at one point, which a raster coder already knows well from the rows above, while removing the cost needs the outline's shape.

The literature agrees. Jumakulyyev and Schultz tried edge-based masks in a PDE-based lossless coder and found that the cost of coding the masks grew disproportionally. In learned image compression the anchor-first checkerboard order costs rate against serial context models, more for weak models than for strong ones (He et al. 2021 and 2022), and multistage orders with few anchors only recover it (Lin et al. 2023). Subscale and multiscale orders lose to raster on small images (Menick and Kalchbrenner 2019, Reed et al. 2017).

Edge detection has two real uses here. The encoder can use it to find objects and send their parameters, and the decoder can adapt its predictor to edges in the pixels it has already decoded.

## What perfect knowledge of the objects is worth to this model

Two oracle probes gave the decoder side information for free.

| archive bytes | seed 3 | seed 2 |
|---|---|---|
| reference | 36,156 | 39,517 |
| additive objects known: the predictor works on the image minus their rendered layer | 35,383 (−773) | 38,699 (−818) |
| one set of predictor weights per true region: background, inside opaque shapes, the rest | 35,837 (−319) | 39,288 (−229) |

Additive objects are all striped discs and about half of the other objects. Knowing them removes 69% of the excess on striped discs and on rings, 28% on lines and 23% on outlines. The background's excess also drops by 88 bytes, because the predictor no longer learns those objects.

Sending objects costs about 8 bytes each for sub-pixel accuracy, about 280 bytes per image, plus a renderer and a fitter in the program, both counted twice. A fitter for six kinds of overlapping objects is unlikely to zip to less than 400 bytes. With this model, object parameters as anchors would roughly break even.

## Candidates, most promising first

Gains are estimated archive bytes per image. Code costs are in bytes of S, which is twice the zipped growth.

### 1. A predictor computed from the generator's spectrum

The background's spectrum is known up to one number, the scale 16/(hi − lo). From the recipe the decoder can compute the covariance for any offset and solve once for the best linear weights, which are the weights kriging uses in geostatistics. Such a predictor needs no learning, does not drift inside objects, and comes with an exact error variance. Added to the mix as a second Gaussian input, it targets the background's whole excess of about 960 bytes, both the cost of learning and the drift caused by objects. The inside of opaque shapes also has a known spectrum, so a second set of weights for it costs a few more bytes. The scale has to be sent in one byte or estimated from the image. Once the variance is computed rather than measured, V − 1/12 becomes the right variance to use.

Estimated gain 300 to 800 bytes for 200 to 300 bytes of S. This is a far cheaper use of the generator than `006_simprior`, which simulated whole images for 960 zipped bytes and saved 119 archive bytes. Here the program needs only the spectrum's formula and a sum of cosines. No published lossless coder seems to do this, since natural images have no known spectrum. The nearest is Gaussian-process regression trained on the image itself (Dai and Xiong 2014).

The same computed weights remove the quadtree's handicap, which is that its coarse levels have few pixels to learn from. On `image.bin` the lab quadtree paid 1.95 bpp for anchors on a step-4 lattice where the ideal is about 1.56, and 1.58 bpp for the next level where the ideal is 1.37. With computed weights every level can reach its ideal, and the calculator says the levels then add up to the raster cost. Anchors would then match raster on the background, and the calculator shows they cannot beat it there.

The first experiment would add the predictor to `lab/costmap.c` and check on clean background that the Gaussian input alone gets within 2% of the ideal error. If it does not, the idea should be dropped.

### 2. Small predictors that forget fast

The reference has one predictor with 112 neighbours and a memory of about 3,300 pixels. Objects are small. An outline crossing is 3 pixels wide, and a striped disc holds a few thousand pixels spread over many rows. A predictor with the 6 to 12 nearest neighbours and a memory of 100 to 300 pixels can lock onto a local edge direction or stripe pattern within one row segment. The image model of paq8px mixes five such predictors with forgetting factors between 0.8 and 0.996. Ulacha and Wernik found that adaptive stages and a bias correction on top of a least-squares predictor are worth about 3% on natural images. Edge-directed prediction (Li and Orchard 2001) gains mostly at edges in the same way.

This targets the outlines, lines, rings and stripes, about 1,750 bytes of excess, and the inside of shapes. Estimated gain 300 to 700 bytes for 150 to 300 bytes of S, since the RLS code can be reused with other sizes. With more inputs the mixer will probably need weights chosen by a small context, because one global weight set has to serve both calm background and edges.

### 3. Keeping objects out of the slow predictor's training

The drift in the background comes from the shared weights learning object interiors. The oracle with one weight set per true region recovered 274 bytes per image. Without side information the decoder has to infer the region from decoded pixels, for example from whether the slow or a fast predictor has done better nearby, and train the slow one only where it wins. Estimated gain 100 to 250 bytes for 40 to 80 bytes of S. Idea 1 removes this problem for the background predictor completely, so this matters only without it.

### 4. Template matching for stripes and lines

A stripe pattern repeats one period back along its direction, and a line repeats one step back along its length. Searching the decoded pixels within about 20 pixels for the best match of a small neighbourhood finds both without knowing the direction. Unno et al. (2019) mixed such a matcher with adaptive predictors and cut the rate of a synthetic shapes image by 27% against MRP. Here it targets about 820 bytes of excess on stripes and lines. Estimated gain 200 to 500 bytes for 250 to 400 bytes of S and about 10 s more CPU per run. This is the least certain of the four.

### 5. Object parameters as anchors

This is the only anchor idea aimed at where the bits are. The encoder finds outlines with an edge detector and fits the generator's six kinds of object to them, circles and lines first with a Hough transform or RANSAC, then a refinement. It sends each object's kind, position, size, angle and level. The decoder renders the anti-aliased masks exactly as `gen_image.py` does. It uses them to subtract additive objects, to choose the predictor for each region, and to predict the outline pixels of opaque objects from both sides.

The oracles put the gain at about 1.1 KB per image with this model, perhaps 1.5 KB with a model built around the masks. Against that come 280 bytes of parameters and the renderer and fitter counted twice. One way to keep the fitter small is to start from the cost map, fit candidates to clusters of expensive pixels, and keep an object only if the archive shrinks by more than its parameters, using the coder itself as the measure. Submitting the compressor and decompressor as separate programs counts the fitter once, but then the model counts three times. The estimated net result lies between −700 and +400 bytes. It is worth trying only after ideas 1 and 2, because a stronger model shrinks what is left for the objects.

### 6. A quadtree refined along edges

This is the idea closest to the original anchor design. Code a coarse lattice, then refine it with interpolation weights estimated from the coarser level in a small window. Li and Orchard's edge-directed interpolation does this, using the coarse level's local covariance to pick the direction. Fine levels see an outline from both sides, and with small stencils hierarchies beat raster. Against it, coarse levels pay for every outline crossing at low resolution, and lines thinner than the lattice step vanish there. Combined with idea 1's computed weights for the background, this is the best chance for a true anchor coder to tie the reference. I expect it to stay a few hundred bytes behind.

## Ideas dropped

- Choosing anchors by edge strength, coding cost or the Laplacian (Belhachmi et al. 2009) and sending their positions. The positions cost more than the anchors can save.
- Orders that follow the data, such as coding the most predictable pixel next. Mentzer et al. (2023) found fixed schedules as good or better, and on the background no order can win.
- Sending low-frequency Fourier coefficients as anchors. The smooth field is already cheap to predict locally, and by the chain rule moving its information to the front saves nothing.
- Correcting the Gaussian's variance for rounding on its own. The mixer already compensates, and the probe saved 1 byte.

## Literature

- Belhachmi, Bucur, Burgeth, Weickert 2009, How to choose interpolation data in images, SIAM J. Appl. Math. https://doi.org/10.1137/080716396
- Mohideen, Peter, Weickert 2021, a comparison of coders for sparse binary inpainting masks. https://arxiv.org/abs/2010.13634
- Peter, Hoffmann, Nedwed, Hoeltgen, Weickert 2016, free masks against tree masks for diffusion-based compression. https://doi.org/10.1016/j.image.2016.05.002
- Jumakulyyev, Schultz 2021, lossless PDE-based compression of 3D medical images, SSVM. https://cg.informatik.uni-bonn.de/backend/v2/files/publications/jumakulyyev-2021-lossless/pdf/jumakulyyev_ssvm21_7db6a35cf7.pdf
- Memon, Neuhoff, Shende 2000, an analysis of some common scanning techniques for lossless image coding, IEEE TIP 9(11). https://doi.org/10.1109/83.877207
- He et al. 2021, checkerboard context model. https://arxiv.org/abs/2103.15306
- He et al. 2022, ELIC. https://openaccess.thecvf.com/content/CVPR2022/papers/He_ELIC_Efficient_Learned_Image_Compression_With_Unevenly_Grouped_Space-Channel_Contextual_CVPR_2022_paper.pdf
- Lin, Sun, Liu, Katto 2023, multistage spatial context models. https://arxiv.org/abs/2302.09263
- Reed et al. 2017, parallel multiscale autoregressive density estimation. https://arxiv.org/abs/1703.03664
- Menick, Kalchbrenner 2019, subscale pixel networks. https://arxiv.org/abs/1812.01608
- Mentzer et al. 2023, M2T. https://arxiv.org/abs/2304.07313
- Dai, Xiong 2014, Gaussian process regression for predictive lossless image coding, Data Compression Conference.
- Li, Orchard 2001, edge-directed prediction for lossless compression of natural images, IEEE TIP 10(6). https://doi.org/10.1109/83.923277
- Li, Orchard 2001, new edge-directed interpolation, IEEE TIP 10(10). https://www.ece.lsu.edu/gunturk/EE7700/NEDI.pdf
- Ulacha, Wernik 2020, lossless image coding with a cascade of predictors, Entropy 22(9). https://www.mdpi.com/1099-4300/22/9/919
- Unno, Kameda, Matsuda, Itoh, Naito 2019, template matching combined with adaptive prediction, EUSIPCO. https://www.eurasip.org/Proceedings/Eusipco/eusipco2019/Proceedings/papers/1570529561.pdf
- Matsuda et al. 2007, MRP, minimum-rate predictors with class maps. https://doi.org/10.1002/scj.20630
- paq8px image model, source code. https://github.com/hxim/paq8px

## Reproduce

```sh
python3 lab/gauss_orders.py 2.5                       # ideal cost of coding orders
mkdir -p /tmp/hr && python3 lab/gen_labeled.py 3 /tmp/hr/s3
cc -O2 -o /tmp/hr/costmap lab/costmap.c -lm
(cd /tmp/hr && ./costmap <s3.bin >a)                  # writes cost.f32, gcost.f32, pred.f32, var.f32
python3 lab/where_bits.py /tmp/hr/s3 /tmp/hr/cost.f32
cp /tmp/hr/s3.tmpl.f32 /tmp/hr/tmpl.f32 && cc -O2 -DTMPL -o /tmp/hr/oracle lab/costmap.c -lm
(cd /tmp/hr && ./oracle <s3.bin | wc -c)              # additive objects known
```
