(ns Guess.core
  (:require (clojure set)
            [clojure.contrib.math :as math])
  (:use [incanter.core :only [mmult matrix factorial view]]
        [incanter.charts :only [function-plot set-y-range add-function]]))

(defn -main
  "I don't do a whole lot."
  [& args]
  (println "Hello, World!"))

;;; ================================================================
;;; Helper functions
;;; ================================================================

(defn fact [n]
  "A factorial that has 0!=1."
  (if (= 0 n)
    1
    (* n (fact (dec n)))))

(defn slope [[x0 y0] [x1 y1]]
  "Given two points [x0 y0] [x1 y1] find the slope."
  (/ (- y0 y1) (- x0 x1)))

(defn line [[x0 y0 :as a] [x1 y1 :as b]]
  "Takes in points a and b where a=[x0 y0] and b=[x1 y1]
   and returns an array of the slope and y-int."
  (let [m (slope a b)]
    [m (- y0 (* m x0))]))

(defn lines [xs]
  "Takes a list of points [p1 p2 p3... ...pn] and returns the
   list of lines ordered as:
        [p1p2 p1p3... ...p1pn p2p3... ...p(n-1)pn]
   where p1p2 is the line passing through p1 and p2"
  (loop [pj    (first xs)
         pks   (rest xs)
         place (rest xs)
         ans []]
    (cond
     (empty? place) ans
     (empty? pks)   (recur (first place) (rest place) (rest place) ans)
     :else          (recur pj (rest pks) place (conj ans (line pj (first pks)))))))


(defn d-eucl [[x0 y0] [x1 y1]]
  "the euclidean distance"
  (math/sqrt (float  (+ (* (- x1 x0) (- x1 x0)) (* (- y1 y0) (- y1 y0))))))

(defn nbhd 
  "Given a metric f (or the Eulidean distance if none given),
    give the three closest points from xs to x."
  ([f x xs]
     (take 3 (sort-by #(apply f [x %]) xs)))
  ([x xs] (nbhd d-eucl x xs)))
    
(defn w-line-dist [[x0 y0] [x1 y1]]
  "Given a line defined by points (xo yo) (x1 y1), give a function
   that gives the distance to parallel line with y-int b."
  (fn [b]
    (* (- (* x1 y0) (* x0 y1) b) (- x1 x0))))

(defn btw [x xs]
  "Given an x and a set of values, gives a vector of the two closest
   points that x falls between, or the two closest points if x is the
   greatest or least value."
  (let [d (split-with #(>= x %) (sort-by #(- % x) xs))]
    (if (empty? (first d))
      (take 2 (last d))
      (if (empty? (last d))
        (reverse (take 2 (reverse (first d))))
        [(last (first d)) (first (last d))]))))

;;; ================================================================
;;; Finding the derivatives and taylor approximations.
;;; ================================================================


(defn deriv-est [a xs]
  "Given a point a and a set of points xs, estimate the
   derivative at that point."
  (let [[b c] (rest (nbhd a xs))
       M (matrix [[1 -1 1] [1 1 -1] [-1 1 1]])
       ans (matrix (map slope [a b a] [b c c]))]
    (first (mmult M ans))))

(defn deriv-ests [xs]
  "Given a set of points, give a list of the points of the
   x-value of each point with the estimated derivatives at that point."
  (map #(apply deriv-est [% xs]) xs))

(defn n-derivs [pts n]
  "Given a set of points return a list of the nth est derivatives."
  (loop [ans pts
         m 0]
    (if (= m n) ans
        (recur (map #(concat (vec %2) [%])
                    (deriv-ests (map #(vector (first %) (last %)) ans)) ans)
               (inc m)))))

(defn taylor [[x0 & fs]]
  "Given a point [x f(x) f'(x)...f^(n)(x)], give the function of the
   Taylor approximation at that point."
  (fn [x]
    (apply + (map-indexed #(/ (* %2 (math/expt (- x x0) %)) (fact %)) fs))))

(defn f-est [pts n]
  "Given a set of points and a number n, returns a set of functions
   that are the estimation of the function at the points given."
  (map taylor (n-derivs pts n)))


;;; ================================================================
;;; Different functions made from the approximations.
;;; ================================================================

(defn f-closest [pts n]
  "Given a set of points, gives a function of x that uses
   the taylor approximation of the closest point to x to the
   nth derivative."
  (fn [x]
    (last (first (sort-by
                  #(math/abs (- (first %) x))
                  (map #(concat % [%2])
                       pts
                       (map #(apply % [x]) (f-est pts n))))))))

(defn f-average [tpts n]
  "Given a set of points, gives a function of x that is the average
   of the taylor approximations to the nth derivative."
  (fn [x]
    (* x (/ (apply + (map #(apply % [x]) (f-est tpts n))) (count tpts)))))

(defn dist-wt [xs x0]
  "Given a set of the x-values of points, and an x0, gives an array of
   weights of each point based on the distance that point is from x0."
  (let [dists  (map #(math/abs (- % x0)) xs)
        D (reduce + dists)]
    (if (some (partial = D) dists)
      (map (fn [x] 1) dists)
      (map #(/ 1 (- 1 (/ % D))) dists))))

(defn f-dist-wt [pts n]
  "Given a set of points pts, gives a function of x that is the weighted
   average of the Taylor approximations to the n-th derivative."
  (fn [x]
    (let [w (dist-wt (map first pts) x)
          f (map #(apply % [x]) (f-est pts n))]
      (apply + (map * w f)))))
  

(defn f-geom-av [pts n]
  "Given a set of points pts, gives a function of x that is the geometric
   average of the Taylor approximations to the n-th derivative."
  (fn [x]
    (math/expt (apply * (map #(apply % [x]) (f-est pts n))) (float (/ 1 (count pts))))))


(defn f-comp [pts n]
  "Given a set of points pts, gives a function of x that is the composition
   of the Taylor approximations to the n-th derivative."
  (fn [x] (reduce #(apply %2 [%]) x (f-est pts n))))


(defn param-line [x [x1 x2]]
  "Given two values x1 and x2, gives a vector [a b] such that
    a*x1 + b*x2 = x , a + b = 1, a,b >= 0."
  (let [t (/ (- x x1) (- x2 x1))]
    [(- 1 t) t]))

(defn f-between [pts n]
  "Given a set of points, gives a function of x that is a weighted
   average of the taylor approximations of the points x lies between
   nth derivative."
  (fn [x]
    (let [[p1 p2 :as p12] (btw x (map first pts))
          [w1 w2 :as ws] (param-line x p12)
          fs (map #(% x) (f-est pts n))]
      (apply +
             (map #(if (= p1 (first %))
                     (* w1 %2)
                     (if (= p2 (first %)) (* w2 %2) 0))
                  pts fs)))))


;;; ================================================================
;;; Stuff to help with charting results.
;;; ================================================================

(defn plotall 
  "Helper function to plot and view a bunch of functions on the same
   chart."
 ([x0 x1 y0 y1 fs]
    (view (set-y-range (reduce #(add-function % %2 x0 x1) 
            (function-plot (first fs) x0 x1)
            (rest fs)) y0 y1)))
 ([x0 x1 fs]
    (view (reduce #(add-function % %2 x0 x1) 
            (function-plot (first fs) x0 x1)
            (rest fs)))))

(defmacro chart [name f method m x0 x1]
  "Helper macro that makes a chart of the guess of f using the method given
   using m random points between x0 and x1 with some random noise."
  `(let [v# (map (fn [x#] (+ (rand ~(- x1 x0)) ~x0)) (repeat ~m 1))
     ~name (~method
                (map (fn [x#] (vector x# ( #(+ (~f %) (rand 1)) x#))) v#))]
     (def ~name (add-function (function-plot ~name ~x0 ~x1) ~f ~x0 ~x1))))


(defn useall []
  "Helper to load all the stuff you need in the Slime-REPL."
  (do (require '[clojure.contrib.math :as math])
                    (use 'Guess.core
                         '[incanter.core :only [mmult matrix factorial view]]
                         '[incanter.charts :only [function-plot
                                                  set-y-range
                                                  add-function]])))
  


    
    






