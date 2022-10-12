package KotlinGeoHash

import com.bahadirarslan.jeodezi.*

fun main() {
  val uh = Coordinate(33.424720, -111.934211)
  val canvas = Coordinate(33.4177786978409, -111.923149300494)
  var greatCircle = GreatCircle()
  var distance = greatCircle.distance(uh, canvas)

  val distanceM = distance * 1000
  println(distanceM)

  // Vincenty test
  var (dist, a1, a2) = vincentyDist(33.424720, -111.934211, 33.4177786978409, -111.923149300494)
  println(dist)
  println(a1)
  println(a2)

  var offset = 1
  while(offset <= 0.1) {
    val newCoord = Coordinate(33.424720, -111.934211 + offset)
    val currDistM = greatCircle.distance(uh, newCoord) * 1000
    val distPerUnit = currDistM / offset
    val oneMeterOffset = 1 / distPerUnit

    println("Offset: $offset, Distance: $currDistM, Unit: $distPerUnit")
    // offset += oneMeterOffset
  }
}