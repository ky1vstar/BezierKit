//
//  Polynomial.swift
//  BezierKit
//
//  Created by Holmes Futrell on 5/15/20.
//  Copyright © 2020 Holmes Futrell. All rights reserved.
//

import CoreGraphics

protocol Polynomial {
    associatedtype Derivative: Polynomial
    func f(_ x: Double, _ scratchPad: UnsafeMutableBufferPointer<Double>) -> Double
    var derivative: Derivative { get }
    var order: Int { get }
    func analyticalRoots(between start: Double, and end: Double) -> [Double]?
}

extension Array: Polynomial where Element == Double {
    typealias Derivative = Self
    var order: Int { return self.count - 1 }
    func split(to x: Double, scratchPad: UnsafeMutableBufferPointer<Double>) -> [Double] {
        _ = f(x, scratchPad)
        let rebased = UnsafeMutableBufferPointer<Double>.init(start: scratchPad.baseAddress!, count: self.count)
        return [Double](rebased)
    }
    func f(_ x: Double, _ scratchPad: UnsafeMutableBufferPointer<Double>) -> Double {
        assert(scratchPad.count >= self.count, "scratchpad will fail here.")
        let count = self.count
        guard count > 0 else { return 0 }
        let oneMinusX = 1.0 - x
        self.withUnsafeBufferPointer { (points: UnsafeBufferPointer<Double>) in
            var i = 0
            while i < count {
                scratchPad[i] = points[i]
                i += 1
            }
            i = 1
            while i < count {
                var j = count-1
                repeat {
                    scratchPad[j] = oneMinusX * scratchPad[j-1] + x * scratchPad[j]
                    j -= 1
                } while j >= i
                i += 1
            }
        }
        return scratchPad[count-1]
    }
    var derivative: [Double] {
        let bufferCapacity = self.order
        guard bufferCapacity > 0 else { return [] }
        let n = Double(bufferCapacity)
        return [Double](unsafeUninitializedCapacity: bufferCapacity) { (buffer: inout UnsafeMutableBufferPointer<Double>, count: inout Int) in
            for i in 0..<bufferCapacity {
                buffer[i] = n * (self[i+1] - self[i])
            }
            count = bufferCapacity
        }
    }
    func analyticalRoots(between start: Double, and end: Double) -> [Double]? {
        let order = self.order
        guard order > 0 else { return [] }
        guard order < 4 else { return nil } // cannot solve
        return Utils.droots(self.map { CGFloat($0) }).compactMap {
            let t = Double($0)
            guard t > start, t < end else { return nil }
            return t
        }
    }
}

private func newton<P: Polynomial>(polynomial: P, derivative: P.Derivative, guess: Double, relaxation: Double = 1, scratchPad: UnsafeMutableBufferPointer<Double>) -> Double {
    let maxIterations = 20
    var x = guess
    for _ in 0..<maxIterations {
        let f = polynomial.f(x, scratchPad)
        guard f != 0.0 else { break }
        let fPrime = derivative.f(x, scratchPad)
        let delta = relaxation * f / fPrime
        let previous = x
        x -= delta
        guard abs(x - previous) > 1.0e-10 else { break }
    }
    return x
}

private func findRootBisection<P: Polynomial>(of polynomial: P, start: Double, end: Double, scratchPad: UnsafeMutableBufferPointer<Double>) -> Double {
    var guess = (start + end) / 2
    var low = start
    var high = end
    let lowSign = polynomial.f(low, scratchPad).sign
    let highSign = polynomial.f(high, scratchPad).sign
    assert(lowSign != highSign)
    let maxIterations = 20
    var iterations = 0
    while high - low > 1.0e-5 {
        let midGuess = (low + high) / 2
        guess = midGuess
        let nextGuessF = polynomial.f(guess, scratchPad)
        if nextGuessF == 0 {
            return guess
        } else if nextGuessF.sign == lowSign {
            low = guess
        } else {
            assert(nextGuessF.sign == highSign)
            high = guess
        }
        iterations += 1
        guard iterations < maxIterations else { break }
    }
    return guess
}

func findRoots(of polynomial: [Double], between start: Double, and end: Double, scratchPad: UnsafeMutableBufferPointer<Double>) -> [Double] {
    assert(start < end)

    var tMin: CGFloat = CGFloat.infinity
    var tMax: CGFloat = -CGFloat.infinity
    var intersected = false

    func x(_ i: Int) -> CGFloat {
        return CGFloat(i) / CGFloat(polynomial.count-1)
    }
    func y(_ i: Int) -> CGFloat {
        return CGFloat(polynomial[i])
    }
    // compute the intersections of each pair of lines with the x axis
    for i in 0..<polynomial.count {
        for j in i+1..<polynomial.count {
            let p1 = CGPoint(x: x(i), y: y(i))
            let p2 = CGPoint(x: x(j), y: y(j))
            let t1 = -p1.y / (p2.y - p1.y)
            guard t1 >= 0, t1 <= 1 else {
                continue
            }
            intersected = true
            let t2 = -p1.y / ((p2.y - p1.y) / (p2.x - p1.x)) + p1.x
            if t2 < tMin {
                tMin = t2
            }
            if t2 > tMax {
                tMax = t2
            }
        }
    }

    guard intersected == true else {
        return [] // no intersections with convex hull
    }

    assert(tMin >= 0 && tMin <= 1)
    assert(tMax >= 0 && tMax <= 1)
    assert(tMax >= tMin)

    // find [adjustedStart, adjustedEnd] range represented by [tMin, tMax] in original polynomial
    func adjustedT(_ t: Double) -> Double {
        return start * (1.0 - t) + end * t
    }
    let adjustedStart = adjustedT(Double(tMin))
    let adjustedEnd = adjustedT(Double(tMax))
    guard adjustedEnd > adjustedStart else {
        return [Double(adjustedStart + adjustedEnd) / 2.0]
    }

    guard tMax - tMin <= 0.8 else {
        // we didn't clip enough of the polynomial off
        // split the polynomial in two and find solutions in each half
        let mid = (start + end) / 2
        let left = polynomial.split(to: 0.5, scratchPad: scratchPad)
        let solutionsLeft = findRoots(of: left, between: start, and: mid, scratchPad: scratchPad)
        let right = [Double](polynomial.reversed()).split(to: 0.5, scratchPad: scratchPad).reversed()
        let solutionsRight = findRoots(of: [Double](right), between: mid, and: end, scratchPad: scratchPad)
        return solutionsLeft + solutionsRight
    }

    // clip the polynomial to [tMin, tMax]
    var clippedPolynomial = polynomial.split(to: Double(tMax), scratchPad: scratchPad)
    let tMinPrime = Double(tMin / tMax)
    clippedPolynomial = [Double]([Double](clippedPolynomial.reversed()).split(to: 1.0 - tMinPrime, scratchPad: scratchPad).reversed())
    return findRoots(of: clippedPolynomial,
                     between: adjustedStart,
                     and: adjustedEnd,
                     scratchPad: scratchPad)
}
