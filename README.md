# Genotyping by sequencing tutorial

authors: Cassandra Elphinstone and Yue Yu \
level: Beginner \
title: GBS tutorial \
date: Jan 2025 

# Background

This is a tutorial to walk you through how to analyze GBS data. We will cover  logging into the Digital Research Alliance servers, demultiplexing, references, mapping, SNP calling and plotting your data.

----

# Programs used in this tutorial



### Perl 

Perl scripts will be used to demultiplex the data

### dDocent 

dDocent is a pipeline that will walk us through building a de novo reference, mapping our reads to that reference and calling SNPs off that reference. 
It uses Rainbow to build the reference,  and FreeBayes or GATK for SNP calling.

### Admixture

Admixture is a program to look at population structure. It is a good way to tell if your data has been demultiplexed correctly. 

### R

R is a programming language that we will use to run a principal component analysis and plot our data.
Libraries needed include: ggplot2, 


-------

# Tutorial style guide

## Text styling

This is how to write _italic text_.

This is how to write **bold text**.

This is how to write **_bold and italic text_**.

Do text superscripts like this 7^th, x^2y or  x^(2y + 3z).


## Lists

### Mixed lists

1. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
2. Integer pharetra arcu ut nisl mollis ultricies.
	* Fusce nec tortor at enim cursus dictum.
	* Phasellus nec urna quis velit eleifend convallis sodales nec augue.
1. In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
1. Nam vitae turpis eu lacus imperdiet mollis id at augue.
1. Sed sed turpis ac dolor mollis accumsan.


## Figures

<figure>
	<a id="fig:example1"></a>
	<img style="width:25%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 1: This figure is 25% of the page width.</figcaption>
</figure>

<figure>
	<a id="fig:example2"></a>
	<img style="width:10%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 2: This figure is only 10% of the page width.</figcaption>
</figure>


# Code

A bit of inline monospaced font can be made `like this`. Larger code blocks can be made by using the code environment:

Java:

```java
public class HelloWorld {

    public static void main(String[] args) {
        // Prints "Hello, World" to the terminal window.
        System.out.println("Hello, World");
    }

}
```

XML:

```xml
	<BirthDeathSkylineModel spec="BirthDeathSkylineModel" id="birthDeath" tree="@tree" contemp="true">
	      <parameter name="origin" id="origin" value ="100" lower="0."/>    
	      <parameter name="R0" id="R0" value="2" lower="0." dimension ="10"/>
	      <parameter name="becomeUninfectiousRate" id="becomeUninfectiousRate" value="1" lower="0." dimension ="10"/>
	      <parameter name="samplingProportion" id="samplingProportion" value="0."/>
	      <parameter name="rho" id="rho" value="1e-6" lower="0." upper="1."/>
	</BirthDeathSkylineModel>
```

R:

```R
	> myString <- "Hello, World!"
	> print (myString)
	[1] "Hello, World!"
```

# Equations

Inline equations: {% eqinline \dot{x} = \sigma(y-x) %}

Displayed equations: 
{% eq \left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right) %}


## Instruction boxes

Use block-quotes for step-by-step instruction that the user should perform (this will produce a framed box on the website):

> The data we have is not the data we want, and the data we need is not the data we have.
> 
> We can input **any** formatted text in here:
>
> - Even
> - Lists
>
> or equations:
>
> {% eq (x_1, \ldots, x_n) \left( \begin{array}{ccc}
      \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
      \vdots & \ddots & \vdots \\
      \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
    \end{array} \right)
  \left( \begin{array}{c}
      y_1 \\
      \vdots \\
      y_n
    \end{array} \right) %}



# Hyperlinks

Add links to figures like this: 

- [Figure 1](#fig:example1) is 25% of the page width.
- [Figure 2](#fig:example2) is 10% of the page width. 

Add links to external URLs like [this](http://www.google.com). 

Links to equations or different sections within the same document are a little buggy.


----

# Useful Links

- Population genomics analyses: [https://github.com/celphin/Population_genomics_Cassiope](https://github.com/celphin/Population_genomics_Cassiope)

----

# Relevant References

{% bibliography --cited --file Tutorial-Template/master-refs.bib %}


