//+build ignore

package main

import (
	. "github.com/mumax/3/engine"
    "math"
    "fmt"
)

func main() {

	defer InitAndClose()()

    cells := 5e-9
    Lx := 1500e-9
    Ly := 700e-9
    Lz := 200e-9

    lx := int(Lx / cells)
    ly := int(Ly / cells)
    lz := int(Lz / cells)
	SetGridSize(lx, ly, lz)
	SetCellSize(cells, cells, cells)

    Msat.Set(4.60545e5)
    Aex.Set(5.72957795e-12)
    Dbulk.Set(0.6e-3)
    // Alpha.Set(0.02)
    LD := 4 * math.Pi * (5.72957795e-12) / (0.6e-3)

    // Magnetisation ---------------------------------------------------------
    
    rbase := Index2Coord(0, 0, 0)
    CENTRE_XSHIFT := Lx * 0.05
    DOT_RAD       := LD / 2
    DOT_LEN       := Ly * 0.75

    // Zi := (float64(lz) * cells) * 0.5
    var rSK [2][2] float64
    rSK[0][0] = -CENTRE_XSHIFT
    rSK[0][1] = 0.0
    rSK[1][0] = CENTRE_XSHIFT
    rSK[1][1] = 0.0
    // X3 := RIGHT_SEP + 5 * SK_XLEN + 2 * SK_SEP

    host := M.Buffer().HostCopy()
    h := host.Vectors()
    n := M.Mesh().Size()

    for iz := 0; iz < n[Z]; iz++ {
        for iy := 0; iy < n[Y]; iy++ {
            for ix := 0; ix < n[X]; ix++ {
                pos := Index2Coord(ix, iy, iz)
                x, y, z := pos[X], pos[Y], pos[Z]

                m := Vector(0.0, -1.0, 0.1)
                for irsk := 0; irsk < 2; irsk++ {

                    xrel := x - rSK[irsk][0]
                    zrel := z - rSK[irsk][1]
                    r_sk := math.Sqrt(xrel * xrel + zrel * zrel)
                    phi_sk := math.Atan2(zrel * 1e9, xrel * 1e9) + 0.5 * math.Pi

                    yrel := math.Abs(rbase[1] + DOT_LEN - y)
                    r_paraboloid := math.Sqrt(yrel * DOT_RAD * DOT_RAD / DOT_LEN)
                    k_sk := math.Pi / r_paraboloid

                    if r_sk < r_paraboloid && y < (rbase[1] + DOT_LEN) {
                        m[0] = math.Sin(k_sk * r_sk) * math.Cos(phi_sk)
                        m[2] = math.Sin(k_sk * r_sk) * math.Sin(phi_sk)
                        m[1] = math.Cos(k_sk * r_sk)
                        break
                    } 
                }
                h[X][iz][iy][ix] = float32(m[X])
                h[Y][iz][iy][ix] = float32(m[Y])
                h[Z][iz][iy][ix] = float32(m[Z])
            }
        }
    }
    M.SetArray(host)

    // -----------------------------------------------------------------------

    SaveAs(&M, "initial_state")

    TableAdd(E_total)
	TableAdd(B_ext)
	TableAdd(E_exch)
	TableAdd(MaxTorque)

    // Relax with conjugate gradient:
    // MinimizerStop = 1e-6
    StopMaxDm = 1e-6
    // Relax with LLG:
    Alpha.Set(0.9)
    Precess = false
    RelaxTorqueThreshold = 1e-3
	// AutoSave(&M, 5e-10)
	TableAutoSave(5e-10)

    Bmax := 210e-3
    Bstep := 20e-3
    simname := fmt.Sprintf("")
    for B := 40e-3; B <= Bmax; B += Bstep {
        // Reset magnetisation:
        M.SetArray(host)
        // simname = fmt.Sprintf("m_initial_By_-%06.0f_mT", B * 1000)
        // SaveAs(&M, simname)

        B_ext.Set(Vector(0, -B, 0))

        // Relax()
        // Relax with specific torque:
        // RunWhile(func() bool { return GetMaxTorque() > 1e-3})
        // Run(12e-9)

        // Relax with LLG and then minimize with SD
        Run(2e-9)
        TableSave()
        Minimize()   // small changes best minimized by minimize()
        TableSave()

        // Save the magnetisation:
        simname = fmt.Sprintf("m_By_-%06.0f_mT", B * 1000)
        SaveAs(&M, simname)
    }
}
