# Changelog

### [0.7.1](https://www.github.com/kjappelbaum/mofchecker/compare/v0.7.0...v0.7.1) (2021-03-26)


### Documentation

* update background on graph hashes ([#117](https://www.github.com/kjappelbaum/mofchecker/issues/117)) ([ec32b54](https://www.github.com/kjappelbaum/mofchecker/commit/ec32b541ba02c809a54546815ae021c241908e0d))

## [0.7.0](https://www.github.com/kjappelbaum/mofchecker/compare/v0.6.0...v0.7.0) (2021-03-24)


### Features

* improve hash ([#108](https://www.github.com/kjappelbaum/mofchecker/issues/108)) ([7765aeb](https://www.github.com/kjappelbaum/mofchecker/commit/7765aeb921ed2d9afd4509a4b4a6687554caf69d))


### Documentation

* fixing links ([c43e919](https://www.github.com/kjappelbaum/mofchecker/commit/c43e919f7095dc0cd182cb8844227eacf6649dce))

## [0.6.0](https://www.github.com/kjappelbaum/mofchecker/compare/v0.5.1...v0.6.0) (2021-03-24)


### Features

* basic zeo++ check ([#104](https://www.github.com/kjappelbaum/mofchecker/issues/104)) ([9bd965a](https://www.github.com/kjappelbaum/mofchecker/commit/9bd965a3d556dadddc44fbe776bdc6521d33b074))


### Bug Fixes

* including requirements.txt in MANIFEST.in ([#109](https://www.github.com/kjappelbaum/mofchecker/issues/109)) ([a960d56](https://www.github.com/kjappelbaum/mofchecker/commit/a960d56dc000c2f4bc70a9f7df4742ea810ca11d))


### Documentation

* adding background to the docs ([5af916b](https://www.github.com/kjappelbaum/mofchecker/commit/5af916b84a7814df6d5d54b756b2469ce6b2ff9f))

### [0.5.1](https://www.github.com/kjappelbaum/mofchecker/compare/v0.5.0...v0.5.1) (2021-03-20)


### Bug Fixes

* dealing with partial occupancy ([#97](https://www.github.com/kjappelbaum/mofchecker/issues/97)) ([98e8357](https://www.github.com/kjappelbaum/mofchecker/commit/98e8357b793a915fc62e4118dc21973cce715ee6))
* graph hash and scaffold hash are now separate properties ([#90](https://www.github.com/kjappelbaum/mofchecker/issues/90)) ([bd76a63](https://www.github.com/kjappelbaum/mofchecker/commit/bd76a6390bc7fce9f7ea367392f413b93f38b94b))
* is_metal choked on D, was not consistent with CSD ([#99](https://www.github.com/kjappelbaum/mofchecker/issues/99)) ([a01519a](https://www.github.com/kjappelbaum/mofchecker/commit/a01519a305cd79fb79dfa915149f56b50bf145ca))

## [0.5.0](https://www.github.com/kjappelbaum/mofchecker/compare/v0.4.0...v0.5.0) (2021-03-17)


### Features

* extend check for undercoordinated metals, closes [#67](https://www.github.com/kjappelbaum/mofchecker/issues/67) ([#70](https://www.github.com/kjappelbaum/mofchecker/issues/70)) ([ceae03b](https://www.github.com/kjappelbaum/mofchecker/commit/ceae03b1afb31d5f02b2c463c836cc2d303365d3))


### Bug Fixes

* fixed graph hash implementation ([#85](https://www.github.com/kjappelbaum/mofchecker/issues/85)) ([46a647f](https://www.github.com/kjappelbaum/mofchecker/commit/46a647f074a9aaf3b7e6ded86ff0e5682faee068))
* typo, closes [#65](https://www.github.com/kjappelbaum/mofchecker/issues/65) ([b43f589](https://www.github.com/kjappelbaum/mofchecker/commit/b43f589ff966a6c8b59a32c8c72a7ec50b630ad7))

## [0.4.0](https://www.github.com/kjappelbaum/mofchecker/compare/v0.3.1...v0.4.0) (2020-12-14)


### Features

* added option to return the Weisfeiler-Lehman hash of the structure graph ([#58](https://www.github.com/kjappelbaum/mofchecker/issues/58)) ([74463bb](https://www.github.com/kjappelbaum/mofchecker/commit/74463bb8f5b4107cf2f0178b92fdf952d8e0ba5d))
* implement check for low metal coordination ([#53](https://www.github.com/kjappelbaum/mofchecker/issues/53)) ([c1676f6](https://www.github.com/kjappelbaum/mofchecker/commit/c1676f632baf54e81d30775897a5fad65b556155))
* implemented charge check ([#52](https://www.github.com/kjappelbaum/mofchecker/issues/52)) ([cd0dbce](https://www.github.com/kjappelbaum/mofchecker/commit/cd0dbce861e206c83081be5d1d4af1f5b7cf7efd))


### Bug Fixes

* patches linear CN=2 case, [#48](https://www.github.com/kjappelbaum/mofchecker/issues/48) ([711cd17](https://www.github.com/kjappelbaum/mofchecker/commit/711cd17eea9951b5acc48aee47f3e83a7807e811))
* underbound N check, closes [#48](https://www.github.com/kjappelbaum/mofchecker/issues/48) ([8c1ba0e](https://www.github.com/kjappelbaum/mofchecker/commit/8c1ba0ef10b9fc6c60b6e155b408cbdb9c722a7c))


### Documentation

* updated features in readme ([3a4be5e](https://www.github.com/kjappelbaum/mofchecker/commit/3a4be5ed2a0caf7b11f388a80029e583ca2bd199))

### [0.3.1](https://www.github.com/kjappelbaum/mofchecker/compare/v0.3.0...v0.3.1) (2020-12-10)


### Documentation

* added installation notes and basic user guide ([4a34d8b](https://www.github.com/kjappelbaum/mofchecker/commit/4a34d8b1b6fd47a3ca385c2c568e0b564dc28cdf))

## [0.3.0](https://www.github.com/kjappelbaum/mofchecker/compare/v0.2.2...v0.3.0) (2020-12-07)


### Features

* expected values and descriptions for checks  ([#42](https://www.github.com/kjappelbaum/mofchecker/issues/42)) ([4f2c33e](https://www.github.com/kjappelbaum/mofchecker/commit/4f2c33e8d36c34c865824df6f3a74b4a5bdc9b08))

### [0.2.2](https://www.github.com/kjappelbaum/mofchecker/compare/v0.2.1...v0.2.2) (2020-12-03)


### Documentation

* Readme links to webapp ([#39](https://www.github.com/kjappelbaum/mofchecker/issues/39)) ([97383be](https://www.github.com/kjappelbaum/mofchecker/commit/97383be40086c64c41b4155aac24ecf8fc56926c))

### [0.2.1](https://www.github.com/kjappelbaum/mofchecker/compare/v0.2.0...v0.2.1) (2020-12-03)


### Documentation

* updated installation section in readme ([443f1b2](https://www.github.com/kjappelbaum/mofchecker/commit/443f1b259f25c75298984f41e23437fd9c1da5d0))
* updated the features in readme ([de653dc](https://www.github.com/kjappelbaum/mofchecker/commit/de653dcf035e193022326c1f50261ca024f3ee87))

## [0.2.0](https://www.github.com/kjappelbaum/mofchecker/compare/v0.1.0...v0.2.0) (2020-12-03)


### Features

* heuristics for missing H on C and N  ([#33](https://www.github.com/kjappelbaum/mofchecker/issues/33)) ([db14cbb](https://www.github.com/kjappelbaum/mofchecker/commit/db14cbb86f69d19cc361b2db20a35e5fe02c17fd))

## 0.1.0 (2020-12-02)


### Features

* added overlap check and property for c/h/overvalent c ([bcf2670](https://www.github.com/kjappelbaum/mofchecker/commit/bcf267048d05b9e5092620f363c148eb0f0cfa92)), closes [#2](https://www.github.com/kjappelbaum/mofchecker/issues/2)
* basic CLI implemented ([d2273fd](https://www.github.com/kjappelbaum/mofchecker/commit/d2273fd6a6cca0bfa15efb242c7b98b7ed281907))
* implemented check for lone atoms and molecules ([fc98d4d](https://www.github.com/kjappelbaum/mofchecker/commit/fc98d4d3de32ffb3ac7910ab2726ebfcc85ad6f5))
* initial commit ([474dab0](https://www.github.com/kjappelbaum/mofchecker/commit/474dab0baa2557fda69a37d1dc8e211ffbebf98f))


### Bug Fixes

* dropping windows support ([5d5e605](https://www.github.com/kjappelbaum/mofchecker/commit/5d5e6055b4098507033c3749c711df0ed878cbee))
* has_carbon in descriptor ([7bba959](https://www.github.com/kjappelbaum/mofchecker/commit/7bba9596df407458875f86bf33b56f3423749c2b))
* installing pytest-cov on github actions. ([ab15d4a](https://www.github.com/kjappelbaum/mofchecker/commit/ab15d4a83b3e4309c98dab1dcb7d24f9adf04bd1))
* subgraph check for newer networkx versions ([85532da](https://www.github.com/kjappelbaum/mofchecker/commit/85532da371eb95ab97fd2551810f0a12478f75cf))


### Documentation

* basic sphinx docs, closes [#27](https://www.github.com/kjappelbaum/mofchecker/issues/27) ([1ae3ae4](https://www.github.com/kjappelbaum/mofchecker/commit/1ae3ae4640b47fa73df4ae582d110dca691d0312))
* basic sphinx docs, closes [#27](https://www.github.com/kjappelbaum/mofchecker/issues/27) ([d5db55e](https://www.github.com/kjappelbaum/mofchecker/commit/d5db55e4b25abce2d0e32c1a4f8cb4e2b2fa5b17))
* updated badges in readme ([40bb8e0](https://www.github.com/kjappelbaum/mofchecker/commit/40bb8e022e0a42cb59953eaa5dc753c52d5b3e7d))
