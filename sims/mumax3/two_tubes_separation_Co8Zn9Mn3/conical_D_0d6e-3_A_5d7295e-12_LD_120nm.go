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

    Ac := 5.72957795e-12
    Dc := 0.6e-3
    Ms := 4.60545e5
    Msat.Set(Ms)
    Aex.Set(Ac)
    Dbulk.Set(Dc)
    // Alpha.Set(0.02)
    LD := 4 * math.Pi * (Ac) / (Dc)
    BD := 0.5 * (Dc * Dc) / (Ac * Ms)

    // Conical phase angle params:
    BDratio := 0.0 / BD

    // Magnetisation ---------------------------------------------------------
    
    // rbase := Index2Coord(0, 0, 0)

    host := M.Buffer().HostCopy()
    h := host.Vectors()
    n := M.Mesh().Size()

    // Set a function to update the magnetisation (using h and host array)
    // using the variables from the main function
    // BDr :: B / BD  ratio used in the Theta angle for the conical background
    magFunction := func(BDr float64) {
        
        for iz := 0; iz < n[Z]; iz++ {
            for iy := 0; iy < n[Y]; iy++ {
                for ix := 0; ix < n[X]; ix++ {
                    pos := Index2Coord(ix, iy, iz)
                    // x, y, z := pos[X], pos[Y], pos[Z]
                    y := pos[Y]

                    m := Vector(0.0, -1.0, 0.1)
                    if math.Abs(BDr) <= 1 {
                        // Set conical BG by default
                        THETAc := math.Acos(BDr)
                        PSIc := 2 * math.Pi * y / LD
                        // fmt.Print(PSIc)
                        m[0] = math.Sin(THETAc) * math.Cos(PSIc)
                        m[1] = math.Cos(THETAc)
                        m[2] = -math.Sin(THETAc) * math.Sin(PSIc)
                    }

                    h[X][iz][iy][ix] = float32(m[X])
                    h[Y][iz][iy][ix] = float32(m[Y])
                    h[Z][iz][iy][ix] = float32(m[Z])
                }
            }
        }
    }

    BDratio = 0.0
    // Specify field ratio
    magFunction(BDratio)
    M.SetArray(host)

    // -----------------------------------------------------------------------


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

    for B := 60e-3; B <= Bmax; B += Bstep {
        
        // Reset magnetisation:
        BDratio = -B / BD
        fmt.Println(BDratio)
        magFunction(BDratio)
        M.SetArray(host)
        // SaveAs(&M, "initial_state")

        // simname = fmt.Sprintf("m_initial_conical_By_-%06.0f_mT", B * 1000)
        // SaveAs(&M, simname)

        B_ext.Set(Vector(0, -B, 0))

        // Relax with LLG and then minimize with SD
        Run(3e-9)
        TableSave()
        Minimize()   // small changes best minimized by minimize()
        TableSave()

        // Save the magnetisation:
        simname = fmt.Sprintf("m_conical_By_-%06.0f_mT", B * 1000)
        SaveAs(&M, simname)
    }
}
