let nint := 1 ;

let {k in int} p1L[k].intervalID := k ;
let {k in int} p1U[k].intervalID := k ;
let {k in int} p2L[k].intervalID := k ;
let {k in int} p2U[k].intervalID := k ;
let {q in int} p1L[q].includeID :=  q;
let {q in int} p1U[q].includeID :=  q;
let {q in int} p2L[q].includeID :=  q;
let {q in int} p2U[q].includeID :=  q;

let p1L[1] := 0.9;
let p1U[1] := 1.1;

#let p1L[2] := 0.9;
#let p1U[2] := 1;

#let p1L[3] := 1;
#let p1U[3] := 1.1;

let p2L[1] := 0.9;
let p2U[1] := 1;
#let p2L[2] := 1;
#let p2U[2] := 1.1;

#let p2L[3] := 1;
#let p2U[3] := 1.1;
#let p1U[2] := 1 ;
#let p1L[3] := 1 ;
#let p1U[3] := 1.05 ;
#let p1L[4] := 1.05 ;
#let p1U[4] := 1.1 ;
#let p1L[5] := 1.05 ;
#let p1U[5] := 1.1 ;
#let p2L[1] := 0.9 ;
#let p2U[1] := 1 ;
#let p2L[2] := 1 ;
#let p2U[2] := 1.1 ;
#let p2L[3] := 1 ;
#let p2U[3] := 1.1 ;
#let p2L[4] := 1 ;
#let p2U[4] := 1.05 ;
#let p2L[5] := 1.05 ;
#let p2U[5] := 1.1 ;