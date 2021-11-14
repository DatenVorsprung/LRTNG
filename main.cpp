/*
 * main.cpp
 *
 ** Lagrangian Reachtubes: The next Generation
 *
 *      Authors: Sophie Gruenbacher, Md Ariful Islam and Jacek Cyranka
 *      Contact: sophie.gruenbacher@tuwien.ac.at
 */

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* IBEX */
#include "ibex.h"

/* Libraries defined by myself */
#include "computeNextStep.h"
#include "computeNorm.h"
#include "ibex_api.h"
#include "processInputFiles.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace ibex;


int main(int argc, char *argv[]) {
    std::time_t timestamp = std::time(0);  // t is an integer type
    bool _DEBUG;
    bool _DEBUG_2 = false;

    /*************** Processing Config and Model related parameter and Dynamics ******************/
    // input: model_file and fdyn_file (e.g. bruss_init.txt and bruss_fdyn.txt)

    double totalTime; // time horizon of computation
    double h; // timestep
    int order; // RK integration order (1, 2 or 4)

    int dim; //dimension of model
    double rad0; //radius of initial box

    bool computeVolume; // compute volume of Reachtube
    bool timeModel; // use model with time variables
    bool output; // print output
    int exactVars;// number of control variables with rad=0

    string model; // files with model parameters and dynamics (e.g. bruss)

    namespace po = boost::program_options;
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("benchmark", po::value<string>(&model)->default_value("bruss"), "choose benchmark")
            ("time_horizon", po::value<double>(&totalTime)->default_value(10.0), "set time horizon")
            ("time_step", po::value<double>(&h)->default_value(0.01), "set timestep")
            ("integration_order", po::value<int>(&order)->default_value(1), "RK integration order (1, 2 or 4)")
            ("compute_volume", po::value<bool>(&computeVolume)->default_value(1), "compute volume of Reachtube (bool)")
            ("print_output", po::value<bool>(&output)->default_value(1), "print output to console (bool)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    cout << "Benchmark model was set to "
         << model << ".\n";
    cout << "Total time horizon was set to "
         << totalTime << ".\n";

    string fdyn_file; // file with c++ function of ODE

    // get parameter values from file
    IntervalVector cx = processInput(totalTime, h, order, dim, rad0, exactVars, computeVolume,
                                     timeModel, output, fdyn_file, model);

    // construct initial radius
    IntervalVector rad(dim, Interval(-rad0, rad0));
    if (exactVars > 0) {
        for (int i = 1; i < exactVars + 1; i++) {
            rad[dim - i] = Interval(0);
        }
    }

    cout << "Initial radius: " << rad[0].ub() << endl << endl;

    IntervalVector rad_new_ellipse(rad);
    IntervalVector rad_circle(rad);

    IntervalVector init_rad(rad);
    IntervalVector init_cx(cx);

    //_DEBUG = output;
    _DEBUG = false;

    //initial xbox
    IntervalVector xbox = cx + rad;
    cout << "Initial Box: " << endl << xbox << endl << endl;

    // get model dynamics from file
    ibex::Function fdyn(fdyn_file.c_str());
    ibex::Function fjacob(fdyn, Function::DIFF);

    if (output)
        std::cout << fdyn << "\n";

    std::vector<ibex::Function> d2fdyn;
    std::vector<std::vector<ibex::Function>> d3fdyn;
    std::vector<std::vector<std::vector<ibex::Function >>> d4fdyn;
    std::vector<std::vector<std::vector<std::vector<ibex::Function >> >> d5fdyn;

    //get symbolic differentiation from fdyn acording to integration order
    symbolicDifferentiation(fdyn, d2fdyn, d3fdyn, d4fdyn, d5fdyn, order, dim);

    /*********** File IO Definition ***************/
    // Output files
    std::ofstream args, ellipse_output, circle_output, usedC_output, usedCi_output,
            xbox_CC_output, xbox_M_output, f_debug, cpp_file, volume_output; // discrete reach-tube
    std::ifstream python_file; //radius output from python file

    string time_stamp = to_string(timestamp);
    string prefix = "saved_outputs/" + time_stamp;

    if (output) {
        ellipse_output.open(prefix + "_ellipse_output.txt",
                            std::ofstream::out); // ellipse: set of time, center, radius and metric
        circle_output.open(prefix + "_circle_output.txt",
                           std::ofstream::out); // circle: set of time, center, radius and metric

        usedC_output.open(prefix + "_usedC_output.txt",
                          std::ofstream::out); // coordinate system switched to in timestep
        usedCi_output.open(prefix + "_usedCi_output.txt",
                           std::ofstream::out); // coordinate system switched to in timestep
        xbox_CC_output.open(prefix + "_xbox_CC_output.txt",
                            std::ofstream::out); // xbox after propagating with interval gradient
        xbox_M_output.open(prefix + "_xbox_M_output.txt",
                           std::ofstream::out); // xbox output in Mi coordinates (usedCxbox)
        volume_output.open(prefix + "_volume_output.txt",
                           std::ofstream::out); // output for initial, current and average volume

        if (_DEBUG) {
            f_debug.open(prefix + "_LRT_debug.txt");
        }
    }

    // set precision of output
    cout.precision(8);

    args.open(prefix + "_args_input.txt", std::ofstream::out); // args with parameters of specific run

    args << "Benchmark model was set to "
         << model << ".\n";
    args << "Total time horizon was set to "
         << totalTime << ".\n";
    args << "Timestep was set to "
         << h << ".\n";

    args.close();

    /****************** Intermediate variable declarations ****************************************/
    int variablesDim = dim;

    if (timeModel) {
        variablesDim -= 1; //dimension of the model excluding the time variable
    }

    ibex::Matrix idMatrix(ibex::Matrix::eye(dim));
    ibex::Matrix idMatrixVariables(ibex::Matrix::eye(variablesDim));

    /* moving forward of gbox with Lohner QR */
    IntervalMatrix gDeltaC(ibex::Matrix::zeros(variablesDim));
    ibex::Matrix gPointC(idMatrixVariables);

    /* Initializing coordinate system with Identity matrix */
    IntervalMatrix usedC(idMatrix);
    IntervalMatrix usedCi(idMatrix); // used coordinate system for the ellipse being formed in this step
    IntervalMatrix oldC(idMatrix);
    IntervalMatrix oldCi(idMatrix); // current coordinate system for the already existing ellipse in this step
    IntervalMatrix initC(idMatrix);
    IntervalMatrix initCi(idMatrix); // coordinate system at time t0 - changes after reinitialization

    /* variables for moving forward of center */
    ibex::Vector cx_error_M0(dim);
    ibex::Vector cx_error_M1(dim);
    ibex::Vector sigma_error_M0(ibex::Vector::zeros(dim));
    ibex::Vector sigma_error_M1(ibex::Vector::zeros(dim));

    /* variables to use Lohner QR */
    IntervalMatrix gbox(idMatrixVariables);
    IntervalMatrix Qg(idMatrixVariables);
    IntervalMatrix Qg_inv(idMatrixVariables);

    /* variables to handle the propagation error of center cbox */
    IntervalMatrix gbox_sigma(variablesDim, variablesDim);
    IntervalVector cbox(cx);
    IntervalMatrix gbox_sigma_C(variablesDim, variablesDim);
    IntervalVector sigma_error_M0_interval(dim);
    IntervalVector sigma_error_M1_interval(dim);

    /* variables of while-loop */
    IntervalMatrix Fmid(variablesDim, variablesDim);
    IntervalMatrix gbox_new_C(variablesDim, variablesDim);
    ibex::Interval cur_sf_circle;
    ibex::Interval cur_sf_ellipse;
    IntervalVector rad_circle_new(dim);
    IntervalVector rad_xbox_M1(dim);
    IntervalVector rad_xbox_M0(dim);
    IntervalVector usedCxbox(dim);
    double rad_max_ellipse;
    double rad_max_circle;

    /* variables for f and fjacob of cx (=cx.mid()), xbox and cbox */
    IntervalMatrix fjacob_cx(dim, dim);
    IntervalMatrix fjacob_xbox(dim, dim);
    IntervalMatrix fjacob_cbox(dim, dim);
    IntervalVector fdyn_cx(dim);
    IntervalVector fdyn_xbox(dim);
    IntervalVector fdyn_cbox(dim);

    /* variables for center stretching factors */
    ibex::Interval cx_sf_circle;
    ibex::Interval cx_sf_ellipse;

    /* variables to compute volumes of gbox and xbox */
    double volumeXbox;
    double sum_volumeXbox = 0;
    double volC = 2. / (double) variablesDim * pow(tgamma((double) variablesDim / 2), -1) *
                  pow(M_PI, ((double) variablesDim / 2)); //volume constant for ellipse
    double volC_ball = pow(tgamma((double) variablesDim / 2 + 1), -1) *
                       pow(M_PI, (double) variablesDim / 2); //volume constant for ball

    /* variables for integration and IST steps */
    double cur_time = 0.0;
    double init_time = cur_time;
    double adapt_h = h; // initialize with user provided value
    int nIter = 1; //start with iteration #1
    int counter_at_init;


    /* Output initial set in the file */

    IntervalMatrix Mmat(dim, dim);

    if (output) {
        ellipse_output << "0 "; // time = 0;
        for (int i = 0; i < variablesDim; i++) {
            ellipse_output << cx[i].mid() << " ";
        }
        ellipse_output << rad[0].ub() << " ";
        IntervalMatrix Mmat = oldC.transpose() * oldC;
        for (int i = 0; i < variablesDim; i++) {
            for (int j = 0; j < variablesDim; j++) {
                ellipse_output << Mmat[i][j].lb() << " ";
            }
        }
        ellipse_output << endl;
    }


    /***** Integration step *****/

    while (cur_time < totalTime) {
        cout << "Integration Iter#" << nIter << endl;

        // final integration step
        if (cur_time + adapt_h > totalTime) {
            adapt_h = totalTime - cur_time;
        }

        cout << "time interval: [" << cur_time << ", " << cur_time + adapt_h << "]" << endl;
        cout << "***************************************************************" << endl << endl << endl;
        if (_DEBUG) {
            f_debug << "time interval: [" << cur_time << ", " << cur_time + adapt_h << "]" << endl;
            f_debug << "*************************************" << endl;
        }
        cout << endl;

        //update current time
        cur_time += adapt_h;


        /* compute once f and fjacob of cx (=cx.mid()), xbox and cbox */
        fjacob_cx = fjacob.eval_matrix(cx);
        fdyn_cx = evalMeanValue(fdyn, fjacob_cx, cx);
        fjacob_xbox = fjacob.eval_matrix(xbox);
        fjacob_cbox = fjacob.eval_matrix(cbox);

        fdyn_xbox = evalMeanValue(fdyn, fjacob_xbox, xbox);
        fdyn_cbox = evalMeanValue(fdyn, fjacob_cbox, cbox);


        /* compute new metric based on last center */
        Fmid = linear_grad(fjacob_cx, adapt_h, dim).submatrix(0, variablesDim - 1, 0, variablesDim -
                                                                                      1); //go between reinit steps & take submatrix if the last variable is time

        ibex::Vector semiAxis(dim);

        getMetric(Fmid, oldC, oldCi, usedC, usedCi, dim, semiAxis, variablesDim);

        /***** compute gradient in coordinate system of used metric to compute stretching factor *****/
        /* version propagating the gradient interval matrix */

        // gbox is a reference, after this function it is already gbox_{t+1}
        // gbox_new_C is already usedC * gbox computed in interval arithmetics and using QR as tight as possible
        // multiplied with initCi to consider also the case where the gradient gets reinitialized
        if (order == 1) {
            gbox_new_C =
                    computeRK1variationalQR(fdyn, fjacob, d2fdyn, fjacob_xbox, fdyn_xbox, xbox, gbox, adapt_h, dim,
                                            gDeltaC, gPointC, Qg, Qg_inv, usedC, idMatrix, variablesDim) *
                    initCi.submatrix(0, variablesDim - 1, 0, variablesDim - 1);
        } else if (order == 2) {
            gbox_new_C = computeRK2variationalQR(fdyn, fjacob, d2fdyn, d3fdyn, fjacob_xbox, fdyn_xbox, xbox, gbox,
                                                 adapt_h, dim,
                                                 gDeltaC, gPointC, Qg, Qg_inv, usedC, idMatrix, variablesDim) *
                         initCi.submatrix(0, variablesDim - 1, 0, variablesDim - 1);
        } else if (order == 4) {
            gbox_new_C =
                    computeRK4variationalQR(fdyn, fjacob, d2fdyn, d3fdyn, d4fdyn, d5fdyn, fjacob_xbox, fdyn_xbox,
                                            xbox, gbox, adapt_h, dim,
                                            gDeltaC, gPointC, Qg, Qg_inv, usedC, idMatrix, variablesDim) *
                    initCi.submatrix(0, variablesDim - 1, 0, variablesDim - 1);
        }


        /* compute stretching factor from time t to t+1 for cbox (cx + propagation error of cx0) */

        // start with gbox_sigma as identity and get next step gbox_sigma in CC and M1 coordinates
        // gbox_sigma_C is the gradient in M1 and gbox_sigma is the one in CC
        gbox_sigma = idMatrixVariables;

        if (order == 1) {
            gbox_sigma_C =
                    computeRK1variationalQR_reinitGbox(fdyn, fjacob, d2fdyn, fjacob_cbox, fdyn_cbox, cbox,
                                                       gbox_sigma,
                                                       adapt_h, dim,
                                                       usedC, idMatrix, variablesDim) *
                    oldCi.submatrix(0, variablesDim - 1, 0, variablesDim - 1);
        } else if (order == 2) {
            gbox_sigma_C =
                    computeRK2variationalQR_reinitGbox(fdyn, fjacob, d2fdyn, d3fdyn, fjacob_cbox, fdyn_cbox, cbox,
                                                       gbox_sigma, adapt_h, dim,
                                                       usedC, idMatrix, variablesDim) *
                    oldCi.submatrix(0, variablesDim - 1, 0, variablesDim - 1);
        } else if (order == 4) {
            gbox_sigma_C =
                    computeRK4variationalQR_reinitGbox(fdyn, fjacob, d2fdyn, d3fdyn, d4fdyn, d5fdyn, fjacob_cbox,
                                                       fdyn_cbox, cbox, gbox_sigma, adapt_h, dim,
                                                       usedC, idMatrix, variablesDim) *
                    oldCi.submatrix(0, variablesDim - 1, 0, variablesDim - 1);
        }


        /***** Compute Stretching factor *****/
        /* using the gradient interval matrix */

        /* SF in CC for the circle */
        cur_sf_circle = get_SF(gbox, variablesDim);

        /* SF in M0,M1 for the ellipse */
        cur_sf_ellipse = get_SF(gbox_new_C, variablesDim);

        /* SF in CC from time t to t+1 for the center box */
        cx_sf_circle = get_SF(gbox_sigma, variablesDim);

        /* SF in M0,M1 from time t to t+1 for the center box */
        cx_sf_ellipse = get_SF(gbox_sigma_C, variablesDim);

        /***** propagate initial center *****/
        cx = RK_center(fdyn, fjacob, fdyn_cx, cx, adapt_h, dim); //cx is already cx.mid() from last step

        if (timeModel)
            cx[dim - 1] = Interval(cur_time);

        /***** compute new radius absorbing error of center propagation *****/
        // error in M0 norm
        cx_error_M0 = cx.rad(); // one-step error of center propagation

        // Corresponding error in M1 norm
        cx_error_M1 = (usedC * (Interval(-1, 1) * cx_error_M0)).mag();

        // absorb current and previous errors to radii
        sigma_error_M0 = cx_sf_circle.ub() * sigma_error_M0 + cx_error_M0;
        sigma_error_M1 = cx_sf_ellipse.ub() * sigma_error_M1 + cx_error_M1;

        // as sigma_error_M0 and sigma_error_M1 describe the same error box just in different
        // coordinate systems, we can minimize them by transforming one to the other and using intersection

        // minimize sigma_error in M0 by taking the minimum with sigma_error_M1 transformed to M0
        sigma_error_M0_interval = Interval(-1, 1) * sigma_error_M0;
        sigma_error_M0_interval &= usedCi * (Interval(-1, 1) *
                                             sigma_error_M1); //intersect with diagonal length of M1-error in cartesian coordinates
        sigma_error_M0 = sigma_error_M0_interval.mag();

        // minimize sigma_error in M1 by taking the minimum with sigma_error_M0 transformed to M1
        sigma_error_M1_interval = Interval(-1, 1) * sigma_error_M1;
        sigma_error_M1_interval &= usedC * (Interval(-1, 1) *
                                            sigma_error_M0); //intersect with diagonal length of M0-error in M1 coordinates
        sigma_error_M1 = sigma_error_M1_interval.mag();

        // add the norm of sigma_error as the radius needed to enclose the box of sigma_error
        rad_circle_new =
                cur_sf_circle.ub() * rad_circle + norm(sigma_error_M0) * IntervalVector(dim, Interval(-1, 1));

        rad_new_ellipse = cur_sf_ellipse.ub() * rad;

        rad_new_ellipse += norm(sigma_error_M1) * IntervalVector(dim, Interval(-1, 1));


        if (timeModel) {
            rad_circle_new[dim - 1] = Interval(0);
            rad_new_ellipse[dim - 1] = Interval(0);
        }

/***** update cbox for the next iteration as tight as possible using intersection *****/
        cbox = cx + Interval(-1, 1) * sigma_error_M0;
        cbox &= cx + usedCi * (Interval(-1, 1) * sigma_error_M1); //intersection

/***** update cx for the next iteration and compute over-approximation xbox of ball in cartesian coordinates *****/

        cx = cx.mid();


/**** Wrap Intersection of ellipse and circle in M0 and M1 boxes as tight as possible ****/

        for (int i = 0; i < dim; i++) {
            rad_xbox_M1[i] =
                    std::min(rad_new_ellipse[i].mag(), rad_circle_new[i].mag() / semiAxis[i]) * Interval(-1,
                                                                                                         1);//rad_new_ellipse is nothing than a factor by which the axis of the ellipse are longer than the columns of usedCi (also in CC),
            // in semiAxis the length of the usedCi-columns are saved
        }

        if (timeModel)
            rad_xbox_M1[dim - 1] = Interval(0);

        usedCxbox = usedC * cx + rad_xbox_M1;

        xbox = usedCi * usedCxbox; // convert xbox into cartesian coordinates

        for (int i = 0; i < dim; i++) {
            rad_xbox_M0[i] =
                    std::min({xbox[i].rad(), rad_circle_new[i].mag(), rad_new_ellipse[i].mag() * semiAxis[i]}) *
                    Interval(-1, 1);// update rad_xbox_M0 intersecting it with the circle in CC,
            // rad_new_ellipse is nothing than a factor by which the axis of the ellipse are longer
            // than the columns of usedCi (also in CC)
        }

        if (timeModel)
            rad_xbox_M0[dim - 1] = Interval(0);

        xbox = xbox.mid() + rad_xbox_M0;


        rad_max_ellipse = rad_new_ellipse.mag().max(); //take the magnitude to be sure that it is the biggest distance from center

        rad_max_circle = rad_circle_new.mag().max();


        // OUTPUTS:
        if (output) {

            /* Write Output set in the file */
            // output of ellipse
            // set represented in (t, cx, r, M)
            // changing to best like in algorithm

            ellipse_output << cur_time << " "; // time > 0;
            for (int i = 0; i < variablesDim; i++) {
                ellipse_output << cx[i].mid() << " ";
            }
            ellipse_output << rad_max_ellipse << " ";
            Mmat = usedC.transpose() * usedC;
            for (int i = 0; i < variablesDim; i++) {
                for (int j = 0; j < variablesDim; j++) {
                    ellipse_output << Mmat[i][j].lb() << " ";
                }
            }
            ellipse_output << endl;

            /* End of Output set in the file */

            /* Write Output set in the file */
            // output of circle

            circle_output << cur_time << " "; // time > 0;
            for (int i = 0; i < variablesDim; i++) {
                circle_output << cx[i].mid() << " ";
            }
            circle_output << rad_max_circle
                          << " ";//take the magnitude because in computation of rad_new it could happen that usedC*cx.mid() is not exactly in the middle
            Mmat = idMatrix;
            for (int i = 0; i < variablesDim; i++) {
                for (int j = 0; j < variablesDim; j++) {
                    circle_output << Mmat[i][j].lb() << " ";
                }
            }
            circle_output << endl;

            /* End of Output set in the file */

            /* Write Output set in the file */
            // output of usedC, which is needed to convert from CC to M1

            usedC_output << cur_time << " "; // time > 0;
            for (int i = 0; i < variablesDim; i++) {
                for (int j = 0; j < variablesDim; j++) {
                    usedC_output << usedC[i][j].lb() << " ";
                }
            }
            usedC_output << endl;

            /* End of Output set in the file */

            /* Write Output set in the file */
            // output of new coordinate System usedCi

            usedCi_output << cur_time << " "; // time > 0;
            for (int i = 0; i < variablesDim; i++) {
                for (int j = 0; j < variablesDim; j++) {
                    usedCi_output << usedCi[i][j].lb() << " ";
                }
            }
            usedCi_output << endl;

            /* End of Output set in the file */


            /* Write Output set in the file */
            // output of ellipse wrapping rectangle in M1 coordinates

            xbox_M_output << cur_time << " "; // time > 0;
            for (int i = 0; i < variablesDim; i++) {
                xbox_M_output << usedCxbox[i].lb() << " ";
                xbox_M_output << usedCxbox[i].ub() << " ";
            }
            xbox_M_output << endl;
            /* End of Output set in the file */


            /* Write Output set in the file */
            // output of ellipse wrapping rectangle in CC coordinates
            xbox_CC_output << cur_time << " "; // time > 0;
            for (int i = 0; i < variablesDim; i++) {
                xbox_CC_output << xbox[i].lb() << " ";
                xbox_CC_output << xbox[i].ub() << " ";
            }
            xbox_CC_output << endl;
            /* End of Output set in the file */
        }


        if (_DEBUG) {
            f_debug << "Center: " << cx << endl;
            f_debug << "SF: " << cur_sf_ellipse.ub() << endl;
            f_debug << "Rad: " << rad_new_ellipse.mag().max() << endl;

            f_debug << "M0: " << endl;
            printIMatrix(oldC.transpose() * oldC, f_debug);

            f_debug << "M1: " << endl;
            printIMatrix(usedC.transpose() * usedC, f_debug);
        }


        if (_DEBUG) {
            f_debug << "xbox: " << xbox << endl;
        }


/**** Updates for the next step ****/

        //update norm for the next iteration
        oldC = usedC;
        oldCi = usedCi;

/** compute sum and average volume of reachtubes **/

        if (computeVolume) {
            double prod_ellipse = 1;
            volumeXbox = 1;

            for (int i = 0; i < variablesDim; i++) {
                prod_ellipse *= rad_max_ellipse * semiAxis[i];
                volumeXbox *= 2 * rad_xbox_M1[i].mag() *
                              semiAxis[i]; // minimum between circle radius and each ellipse semi-axis
            }

            double volumeEllipse = volC * prod_ellipse;
            double volumeBall;
            volumeBall = volC_ball * pow(rad_max_circle, variablesDim);


            if (_DEBUG_2) {
                cout << "volumeEllipse: " << volumeEllipse << endl << endl;
                    cout << "volumeBall: " << volumeBall << endl << endl;
                cout << "volumeXbox: " << volumeXbox << endl << endl;
            }

            volumeXbox = min(volumeEllipse, volumeXbox);
            volumeXbox = min(volumeBall, volumeXbox);
            sum_volumeXbox += volumeXbox;

            if (output) {
                cout << endl << "Volume:	" << volumeXbox << endl << endl;
                cout << endl << "Average volume:	" << sum_volumeXbox / (nIter) << endl << endl;

                /* Write Output set in the file */
                // output of volume
                volume_output << "Iterations: " << nIter << "; ";
                volume_output << "Current time: " << cur_time << "; ";
                volume_output << "Average volume: " << (sum_volumeXbox / (nIter)) << "; ";
                volume_output << "Current Volume: " << volumeXbox << "; ";
                volume_output << endl;
                /* End of Output set in the file */
            }
        }

        nIter++;
    }// END INTEGRATION LOOP
    nIter--;

    if (computeVolume) {
        double initialVolume = volC_ball * pow(rad0, variablesDim);
        cout << endl << endl << endl << endl << "Volume after " << nIter << " iterations:" << endl;
        cout << endl << "Average Volume:	" << (sum_volumeXbox / (nIter)) << endl;
        cout << endl << "Initial Volume:	" << initialVolume << endl;
        cout << endl << "Final Volume:	" << volumeXbox << endl;
    }

    cout << endl << "Plot result by running:" << endl;
    cout << endl << "cd plot_results" << endl;
    cout << "source venv/bin/activate" << endl;
    cout << "python plot_LRTNG.py";
    cout << " --dim " << variablesDim << " --time_horizon " << cur_time << " --time_step " << h;
    cout << " --circle_file " << time_stamp << "_circle_output.txt";
    cout << " --ellipse_file " << time_stamp << "_ellipse_output.txt" << endl;
    cout << "deactivate" << endl;
    cout << "cd .." << endl << endl << endl;

    /* close output file streams */
    if (output) {
        xbox_CC_output.close();
        xbox_M_output.close();
        ellipse_output.close();
        circle_output.close();
        usedC_output.close();
        usedCi_output.close();
        volume_output.close();
        if (_DEBUG)
            f_debug.close();
    }

    return EXIT_SUCCESS;
}