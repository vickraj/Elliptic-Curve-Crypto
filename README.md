# Elliptic-Curve-Crypto
A suite of Algorithms for computations on Elliptic Curves over prime fields.
Currently this supports:
*ECDSA 
*Primality testing 
*Factoring 
*Discrete Logs on elliptic curves
*Weil pairings 
*One round Tripartite Diffie Hellman Key Exchange
*Arithmetic on elliptic curves (including random point generation, scaling, adding, etc.)
*Prime field arithmetic (including inverses, square roots (Tonelli-Shanks))

## Overview:
This suite was built as a way of testing my understanding of several algorithms that I learned the theory for. I was interested in seeing the details of each algorithm in action, in a way where I could play around with several assumptions or try my own values or modify them to my hearts desire. Also, since I wanted to implement the algorithms myself, I refrained from using Sage (which honestly made this both way more annoying and enlightening). 

## Dependencies:
sympy -- calculating the weil pairing in particular, among other things, uses sympy. 
numpy -- several calculations use numpy's pseudo random number generator. 

## Usage:
The ecc.py file has examples of all of the other modules being used, as well as a small terminal that allows one to use these.

## TODO:
*Include support for prime power fields. Will likely learn a bit more about Programming Languages and create my own types for this. 

*Include Schoof's algorithm or another version of efficient point counting, which is necessary to do ECPP and the MOV attack.

*Include support for larger integers/primes.

*Other TODO's are in each of the modules separately. 
