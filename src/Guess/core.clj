(ns Guess.core)

(defn -main
  "I don't do a whole lot."
  [& args]
  (println "Hello, World!"))

(defn slope [a b]
  (let [x0 (a 0) y0 (a 1)
        x1 (b 0) y1 (b 1)]
    (/ (- y0 y1) (- x0 x1))))

(defn line [a b]
  (let [m (slope a b)
        x0 (a 0)
        y0 (a 1)]
    [m (+ (* m x0 -1) y0)]))


(defn lines [xs]
  (loop [pj    (first xs)
         pks   (rest xs)
         place (rest xs)
         ans []]
    (cond
     (empty? place) ans
     (empty? pks)   (recur (first place) (rest place) (rest place) ans)
     :else          (recur pj (rest pks) place (conj ans (line pj (first pks)))))))


(defn graph [f xs]
  (loop [ps xs ans []]
    (if (empty? ps)
      ans
      (recur (rest ps) (conj ans (f (first ps)))))))

(defn solve [l1 l2]
  (let [m1 (l1 0) b1 (l1 1)
        m2 (l2 0) b2 (l2 1)
        x (/ (- b2 b1) (- m1 m2))
        y1 (+ (* x m1) b1)
        y2 (+ (* x m2) b2)]
    (if (= y1 y2)
      [x y1]
      :error)))