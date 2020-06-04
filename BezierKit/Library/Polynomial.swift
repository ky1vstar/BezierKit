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
    func derivative(_: UnsafeMutableBufferPointer<Double>) -> Derivative
    var order: Int { get }
    func analyticalRoots(between start: Double, and end: Double, scratchPad: UnsafeMutableBufferPointer<Double>) -> UnsafeMutableBufferPointer<Double>?
}

internal extension Polynomial {
    func analyticalRoots(between start: Double, and end: Double) -> [Double]? {
        let scratchPad = UnsafeMutableBufferPointer<Double>.allocate(capacity: 3)
        defer { scratchPad.deallocate() }
        guard let result = self.analyticalRoots(between: start, and: end, scratchPad: scratchPad) else {
            return nil
        }
        return [Double](result)
    }
}

extension UnsafeMutableBufferPointer {
    mutating func push(_ count: Int) -> UnsafeMutableBufferPointer<Element> {
        assert(self.count >= count)
        let value = UnsafeMutableBufferPointer(start: self.baseAddress, count: count)
        self = UnsafeMutableBufferPointer<Element>(start: self.baseAddress! + count, count: self.count - count)
        return value
    }
    mutating func pop(_ count: Int) {
        self = UnsafeMutableBufferPointer<Element>(start: self.baseAddress! - count, count: self.count + count)
    }
}

extension UnsafeMutableBufferPointer: Polynomial where Element == Double {
    typealias Derivative = Self
    var order: Int { return self.count - 1 }
    func f(_ x: Double, _ scratchPad: UnsafeMutableBufferPointer<Double>) -> Double {
        let count = self.count
        guard count > 0 else { return 0 }
        let oneMinusX = 1.0 - x
        var i = 0
        while i < count {
            scratchPad[i] = self[i]
            i += 1
        }
        i = count - 1
        while i > 0 {
            var j = 0
            repeat {
                scratchPad[j] = oneMinusX * scratchPad[j] + x * scratchPad[j+1]
                j += 1
            } while j < i
            i -= 1
        }
        return scratchPad[0]
    }
    func derivative(_ buffer: UnsafeMutableBufferPointer<Double>) -> UnsafeMutableBufferPointer<Double> {
        let bufferCapacity = self.order
        assert(buffer.count >= bufferCapacity)
        assert(bufferCapacity > 0)
        let n = Double(bufferCapacity)
        for i in 0..<bufferCapacity {
            buffer[i] = n * (self[i+1] - self[i])
        }
        return UnsafeMutableBufferPointer(start: buffer.baseAddress!, count: bufferCapacity)
    }
    func analyticalRoots(between start: Double, and end: Double, scratchPad: UnsafeMutableBufferPointer<Double>) -> UnsafeMutableBufferPointer<Double>? {
        var scratchPad = scratchPad
        let order = self.order
        guard order > 0 else { return scratchPad.push(0) }
        guard order < 4 else { return nil } // cannot solve
        var i = 0
        func c(_ x: CGFloat) {
            scratchPad[i] = Double(x)
            i += 1
        }
        switch self.count {
        case 4:
            Utils.droots(CGFloat(self[0]), CGFloat(self[1]), CGFloat(self[2]), CGFloat(self[3]), callback: c)
        case 3:
            Utils.droots(CGFloat(self[0]), CGFloat(self[1]), CGFloat(self[2]), callback: c)
        case 2:
            Utils.droots(CGFloat(self[0]), CGFloat(self[1]), callback: c)
        default:
            assertionFailure("?")
        }
        scratchPad = UnsafeMutableBufferPointer(start: scratchPad.baseAddress!, count: i)
        return scratchPad
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

func findRoots(of polynomial: UnsafeMutableBufferPointer<Double>, between start: Double, and end: Double, roots: inout UnsafeMutableBufferPointer<Double>, scratchPad: UnsafeMutableBufferPointer<Double>) {
    assert(start < end)

    var scratchPad = scratchPad

    let analyticalRoots: UnsafeMutableBufferPointer<Double> = scratchPad.push(3)
    defer { scratchPad.pop(3) }

    if let analyticalRoots = polynomial.analyticalRoots(between: start, and: end, scratchPad: scratchPad) {
        roots = UnsafeMutableBufferPointer<Double>(start: roots.baseAddress, count: analyticalRoots.count)
        for i in 0..<analyticalRoots.count {
            roots[i] = analyticalRoots[i]
        }
        return
    }

    var derivativeStorage = scratchPad.push(polynomial.count - 1)
    defer { scratchPad.pop(derivativeStorage.count) }
    let derivative = polynomial.derivative(derivativeStorage)

    var criticalPoints = scratchPad.push(derivative.order)
    defer { scratchPad.pop(derivative.order) }
    findRoots(of: derivative, between: start, and: end, roots: &criticalPoints, scratchPad: scratchPad)

    let intervals = scratchPad.push(criticalPoints.count + 2)
    defer { scratchPad.pop(criticalPoints.count + 2) }
    intervals[0] = start
    intervals[intervals.count - 1] = end
    for i in 0..<criticalPoints.count {
        intervals[i+1] = criticalPoints[i]
    }

    var lastFoundRoot: Double?

    var rootsFound = 0

    for i in (0..<intervals.count-1) {
        let start   = intervals[i]
        let end     = intervals[i+1]
        let fStart  = polynomial.f(start, scratchPad)
        let fEnd    = polynomial.f(end, scratchPad)
        let root: Double
        if fStart * fEnd < 0 {
            // TODO: if a critical point is a root we take this
            // codepath due to roundoff and  converge only linearly to one end of interval
            let guess = (start + end) / 2
            let newtonRoot = newton(polynomial: polynomial, derivative: derivative, guess: guess, scratchPad: scratchPad)
            if start < newtonRoot, newtonRoot < end {
                root = newtonRoot
            } else {
                // newton's method failed / converged to the wrong root!
                // rare, but can happen roughly 5% of the time
                // see unit test: `testDegree4RealWorldIssue`
                root = findRootBisection(of: polynomial, start: start, end: end, scratchPad: scratchPad)
            }
        } else {
            let guess = end
            let value = newton(polynomial: polynomial, derivative: derivative, guess: guess, scratchPad: scratchPad)
            guard abs(value - guess) < 1.0e-5 else {
                continue // did not converge near guess
            }
            guard abs(polynomial.f(value, scratchPad)) < 1.0e-10 else {
                continue // not actually a root
            }
            root = value
        }
        if let lastFoundRoot = lastFoundRoot {
            guard lastFoundRoot + 1.0e-5 < root else {
                continue // ensures roots are unique and ordered
            }
        }
        lastFoundRoot = root
        roots[rootsFound] = root
        rootsFound += 1
    }
    roots = UnsafeMutableBufferPointer<Double>(start: roots.baseAddress!, count: rootsFound)
}
