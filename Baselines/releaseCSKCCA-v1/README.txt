% %
Created by: Zhiqiang Li
Email: lzq115@163.com


The code for the paper of
Zhiqiang Li, Xiao-Yuan Jing, Fei Wu, Xiaoke Zhu, Baowen Xu and Shi Ying. 
"Cost-Sensitive Transfer Kernel Canonical Correlation Analysis for Heterogeneous Defect Prediction" 
which has been published to the Joural of Automated Software Engineering.

Please kindly cite our paper if you would like to use the code for your research.

%%%
Running demo_CSTKCCA

MATLAB R2014a, 64bit operating system (the liblinear just only provides ".mexw64" files, you can remake these files according to the README file)
%%%


Note, we use km_kernel and km_kernel subfunctions, which are part of the Kernel Methods Toolbox
of Steven Van Vaerenbergh.
@book{vaerenbergh2010kernel,
	title={Kernel methods for nonlinear identification, equalization and separation of signals},
	author={Vaerenbergh, Steven van},
	year={2010},
	publisher={Universidad de Cantabria}
}

And we use LR classifier from LIBLINEAR, which is a library for
large-scale regularized linear classification and regression
(http://www.csie.ntu.edu.tw/~cjlin/liblinear). It is very easy to use
as the usage and the way of specifying parameters are the same as that
of LIBLINEAR.
@article{fan2008liblinear,
	title={LIBLINEAR: A library for large linear classification},
	author={Fan, Rong-En and Chang, Kai-Wei and Hsieh, Cho-Jui and Wang, Xiang-Rui and Lin, Chih-Jen},
	journal={The Journal of Machine Learning Research},
	volume={9},
	pages={1871--1874},
	year={2008},
	publisher={JMLR. org}
}

%%% NOTE %%%
The software is free for academic use only, and shall not be used, rewritten, or adapted as the basis of a commercial product without first obtaining permission from the authors. The authors make no representations about the suitability of this software for any purpose. It is provided "as is" without express or implied warranty.

