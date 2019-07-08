# A Quantum program to illustrate the ISING model which is NP-COMPLETE for 72-Qubit Quantum Computer ![](https://travis-ci.com/CodeBreaker444/A-Quantum-program-to-illustrate-the-ISING-model-which-is-NP-COMPLETE.svg?branch=master)
> `Note`: Before anything *Cirq* which is a quantum library for python (A python framework for creating, editing, and invoking Noisy Intermediate Scale Quantum (NISQ) circuits.) is still in alpha stage and is highly susceptible to syntax changes or even basic framework changes. So, only use the version specified in the `requirements.txt`.

> `Note 2`: As travis builds halts if execution time exceeds 10 Mins. I have reduced the circuit_repeats to 1 but it's an obligatory value to be at 100 repetitions. Please, change it if used for any official research purpose.

## Requirements:
- Python3
- pycharm(Optional)

> Note: This is a small quantum computer program just to prove a theory and at the moment writing this readme their is no quantum processor which is capable of doing classical or large quantum computations.
## Installation
``` pip3 install -r requirements.txt ```
## Packages Used
```
import cirq
import random
import sympy
import numpy as np
import timeit

```

## You will see something like
```
Moment 0: H((0, 0)) and H((0, 2)) and H((1, 1)) and H((2, 0)) and H((2, 2))
Moment 1: X((0, 1)) and X((1, 0)) and X((1, 2)) and X((2, 1))
                           ┌──────┐   ┌──────────┐       ┌──────────┐
(0, 0): ───X^0.1───X───────────────────@─────────────X────X─────────────@───────X───────────────────────────────
                                       │                                │
(0, 1): ───X^0.1───Z^0.2────X──────────┼────@────────X────X─────────────@^0.3───X───────@───────────────────────
                                       │    │                                           │
(0, 2): ───X^0.1───Z^0.2─────@─────────┼────┼───────────────────────────────────────────@^0.3───────────────────
                             │         │    │
(1, 0): ───X^0.1───Z^0.2────X┼─────────@^0.3┼────────X────@─────────────@───────────────────────────────────────
                             │              │             │             │
(1, 1): ───X^0.1───Z^0.2────X┼──────────────@^0.3────X────┼────@────────@^0.3───X───────@───────X───────────────
                             │                            │    │                        │
(1, 2): ───X^0.1─────────────@^0.3─────@─────────────X────┼────┼────────────────────────@^0.3───X───────────────
                                       │                  │    │
(2, 0): ───X^0.1───────────────────────┼──────────────────@^0.3┼────────X───────@───────X───────────────────────
                                       │                       │                │
(2, 1): ───X^0.1───Z^0.2───────────────┼───────────────────────@^0.3────X───────@^0.3───X───────X───@───────X───
                                       │                                                            │
(2, 2): ───X^0.1───────────────────────@^0.3─────────X──────────────────────────────────────────────@^0.3───X───
                           └──────┘   └──────────┘       └──────────┘
transverse fields: [[-1, 1, 1], [1, -1, -1], [-1, 1, -1]]
row j fields: [[1, 1, 1], [1, 1, 1]]
column j fields: [[1, 1], [1, -1], [-1, -1]]
(0, 0): ───X^0.1───

(0, 1): ───X^0.1───

(1, 0): ───X^0.1───

(1, 1): ───X^0.1───
Counter({0: 1})
Counter({5: 1})
Value of the objective function 5.0
                           ┌───────────────┐   ┌───────────────┐                                                                       ┌─────────────────────┐   ┌─────────────────────┐
(0, 0): ───X^0.1────────────@───────────────────@──────────────────────────────────────────────────────────M('x')───X^alpha─────────────@─────────────────────────@────────────────────────────────────────────────────────────────────────M('x')───
                            │                   │                                                          │                            │                         │                                                                        │
(0, 1): ───X^0.1───Z^0.2────┼────@──────────────@^0.3──────────────@───────────────────────────────────────M────────X^alpha───Z^beta────┼──────@──────────────────@^gamma──────────────────@───────────────────────────────────────────────M────────
                            │    │                                 │                                       │                            │      │                                           │                                               │
(0, 2): ───X^0.1───Z^0.2────┼────┼────@────────────────────────────@^0.3───────────────────────────────────M────────X^alpha───Z^beta────┼──────┼──────@────────────────────────────────────@^gamma─────────────────────────────────────────M────────
                            │    │    │                                                                    │                            │      │      │                                                                                    │
(1, 0): ───X^0.1───Z^0.2────@^0.3┼────┼─────────@──────────────────@───────────────────────────────────────M────────X^alpha───Z^beta────@^gamma┼──────┼───────────@────────────────────────@───────────────────────────────────────────────M────────
                                 │    │         │                  │                                       │                                   │      │           │                        │                                               │
(1, 1): ───X^0.1─────────────────@^0.3┼─────────┼────@─────────────@^0.3───X───────@───────X───────────────M────────X^alpha────────────────────@^gamma┼───────────┼──────@─────────────────@^gamma───X─────────@─────────X─────────────────M────────
                                      │         │    │                             │                       │                                          │           │      │                                     │                           │
(1, 2): ───X^0.1──────────────────────@^0.3─────┼────┼────@────────X───────────────@^0.3───X───────────────M────────X^alpha───────────────────────────@^gamma─────┼──────┼──────@──────────X───────────────────@^gamma───X─────────────────M────────
                                                │    │    │                                                │                                                      │      │      │                                                          │
(2, 0): ───X^0.1────────────────────────────────@^0.3┼────┼────────X───────@───────X───────────────────────M────────X^alpha───────────────────────────────────────@^gamma┼──────┼──────────X─────────@─────────X───────────────────────────M────────
                                                     │    │                │                               │                                                             │      │                    │                                     │
(2, 1): ───X^0.1───Z^0.2─────────────────────────────@^0.3┼────────X───────@^0.3───X───────X───@───────X───M────────X^alpha───Z^beta─────────────────────────────────────@^gamma┼──────────X─────────@^gamma───X─────────X───@─────────X───M────────
                                                          │                                    │           │                                                                    │                                            │             │
(2, 2): ───X^0.1──────────────────────────────────────────@^0.3────X───────────────────────────@^0.3───X───M────────X^alpha─────────────────────────────────────────────────────@^gamma────X─────────────────────────────────@^gamma───X───M────────
                           └───────────────┘   └───────────────┘                                                                       └─────────────────────┘   └─────────────────────┘
Minimum objective value is -12.0.
Time:  26.57741808299761

```
> `Min. Sweep value changes due to unpredictable quantum errors but ISING model problem can be solved`
> `Check it in Action:` [Click Here](https://travis-ci.com/CodeBreaker444/A-Quantum-program-to-illustrate-the-ISING-model-which-is-NP-COMPLETE/builds/118221709)
## Comments are present at important stages for better view of its working which avoids blind execution :)
## Personal INFO:
`Donations Help Me to Keep The Support and Development:` [Click Here](https://paypal.me/zer0error).

`FollowMe:` [Click Here](https://facebook.com/zer0error/).

`Google Play:` [Codebreaker](https://play.google.com/store/apps/dev?id=8331274631553271784&hl=en).

`Website:` [Personal](https://govardhanchitrada.me).
