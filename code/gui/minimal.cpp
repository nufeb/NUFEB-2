/////////////////////////////////////////////////////////////////////////////
// Name:        minimal.cpp
// Purpose:     Minimal wxWidgets sample
// Author:      Julian Smart
// Modified by:
// Created:     04/01/98
// RCS-ID:      $Id$
// Copyright:   (c) Julian Smart
// Licence:     wxWindows licence
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWidgets headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include <wx/grid.h>

// ----------------------------------------------------------------------------
// resources
// ----------------------------------------------------------------------------

// the application icon (under Windows and OS/2 it is in resources and even
// though we could still include the XPM here it would be unused)
#ifndef wxHAS_IMAGES_IN_RESOURCES
    #include "sample.xpm"
#endif

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

// Define a new application type, each program should derive a class from wxApp
class MyApp : public wxApp
{
public:
    // override base class virtuals
    // ----------------------------

    // this one is called on application startup and is a good place for the app
    // initialization (doing it here and not in the ctor allows to have an error
    // return: if OnInit() returns false, the application terminates)
    virtual bool OnInit();
};

// Define a new frame type: this is going to be our main frame
class MyFrame : public wxFrame
{
public:
    // ctor(s)
    MyFrame(const wxString& title);

    // event handlers (these functions should _not_ be virtual)
    void OnQuit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
    void Browse(wxCommandEvent& event);

private:
    // any class wishing to process wxWidgets events must use this macro
    wxDECLARE_EVENT_TABLE();

    wxTextCtrl *inputFile;
};

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

// IDs for the controls and the menu commands
enum
{
    // menu items
    Minimal_Quit = wxID_EXIT,

    // it is important for the id corresponding to the "About" command to have
    // this standard value as otherwise it won't be handled properly under Mac
    // (where it is special and put into the "Apple" menu)
    Minimal_About = wxID_ABOUT,
    BUTTON_Browse = wxID_HIGHEST + 1
};

// ----------------------------------------------------------------------------
// event tables and other macros for wxWidgets
// ----------------------------------------------------------------------------

// the event tables connect the wxWidgets events with the functions (event
// handlers) which process them. It can be also done at run-time, but for the
// simple menu events like this the static method is much simpler.
wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(Minimal_Quit,  MyFrame::OnQuit)
    EVT_MENU(Minimal_About, MyFrame::OnAbout)
    EVT_BUTTON ( BUTTON_Browse, MyFrame::Browse )
wxEND_EVENT_TABLE()

// Create a new application object: this macro will allow wxWidgets to create
// the application object during program execution (it's better than using a
// static object for many reasons) and also implements the accessor function
// wxGetApp() which will return the reference of the right type (i.e. MyApp and
// not wxApp)
IMPLEMENT_APP(MyApp)

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// the application class
// ----------------------------------------------------------------------------

// 'Main program' equivalent: the program execution "starts" here
bool MyApp::OnInit()
{
    // call the base class initialization method, currently it only parses a
    // few common command-line options but it could be do more in the future
    if ( !wxApp::OnInit() )
        return false;

    // create the main application window
    MyFrame *frame = new MyFrame("NUFEB Model");

    // and show it (the frames, unlike simple controls, are not shown when
    // created initially)
    frame->Show(true);

    // success: wxApp::OnRun() will be called which will enter the main message
    // loop and the application will run. If we returned false here, the
    // application would exit immediately.
    return true;
}

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

// frame constructor
MyFrame::MyFrame(const wxString& title)
       : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(2000,1000))
{
    // set the frame icon
    SetIcon(wxICON(sample));

#if wxUSE_MENUS
    // create a menu bar
    wxMenu *fileMenu = new wxMenu;

    // the "About" item should be in the help menu
    wxMenu *helpMenu = new wxMenu;
    helpMenu->Append(Minimal_About, "&About\tF1", "Show about dialog");

    fileMenu->Append(Minimal_Quit, "E&xit\tAlt-X", "Quit this program");

    // now append the freshly created menu to the menu bar...
    wxMenuBar *menuBar = new wxMenuBar();
    menuBar->Append(fileMenu, "&File");
    menuBar->Append(helpMenu, "&Help");

    // ... and attach this menu bar to the frame
    SetMenuBar(menuBar);
#endif // wxUSE_MENUS

#if wxUSE_STATUSBAR
    // create a status bar just for fun (by default with 1 pane only)
    CreateStatusBar(2);
    SetStatusText("Enter parameters");
#endif // wxUSE_STATUSBAR


    wxPanel *panel = new wxPanel(this, -1);

    wxBoxSizer *main = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer *sizer = new wxBoxSizer(wxVERTICAL);


    sizer->Add(new wxStaticText(
        panel, -1, "All units should be specified in SI"));

    sizer->Add(new wxStaticText(
        panel, -1, "", wxDefaultPosition));

    wxBoxSizer *vbox1 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox2 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox3 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox4 = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer *hboxBoundary = new wxBoxSizer(wxHORIZONTAL);


    wxRadioButton *periodicX = new wxRadioButton(panel, -1, 
        "", wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    wxRadioButton *fixedX = new wxRadioButton(panel, -1, 
        "");
    wxRadioButton *periodicY = new wxRadioButton(panel, -1, 
        "", wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    wxRadioButton *fixedY = new wxRadioButton(panel, -1, 
        "");
    wxRadioButton *periodicZ = new wxRadioButton(panel, -1, 
        "", wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    wxRadioButton *fixedZ = new wxRadioButton(panel, -1, 
        "");

    vbox1->Add(new wxStaticText(
        panel, -1, "Boundary"));

    vbox1->Add(new wxStaticText(
        panel, -1, "Periodic"));

    vbox1->Add(new wxStaticText(
        panel, -1, "Fixed"));

    vbox2->Add(new wxStaticText(
        panel, -1, "x"), 0, wxALIGN_CENTER);

    vbox2->Add(periodicX, 0, wxALIGN_CENTER);

    vbox2->Add(fixedX, 0, wxALIGN_CENTER);

    vbox3->Add(new wxStaticText(
        panel, -1, "y"), 0, wxALIGN_CENTER);

    vbox3->Add(periodicY, 0, wxALIGN_CENTER);

    vbox3->Add(fixedY, 0, wxALIGN_CENTER);

    vbox4->Add(new wxStaticText(
        panel, -1, "z"), 0, wxALIGN_CENTER);

    vbox4->Add(periodicZ, 0, wxALIGN_CENTER);

    vbox4->Add(fixedZ, 0, wxALIGN_CENTER);

    hboxBoundary->Add(vbox1);
    hboxBoundary->Add(vbox2);
    hboxBoundary->Add(vbox3);
    hboxBoundary->Add(vbox4);

    sizer->Add(hboxBoundary);

    sizer->Add(new wxStaticText(
        panel, -1, "", wxDefaultPosition));


    wxBoxSizer *hboxInput = new wxBoxSizer(wxHORIZONTAL);


    inputFile = new wxTextCtrl(panel, -1, "", wxDefaultPosition, wxSize(400,-1),  
    wxTE_LEFT | wxTE_READONLY);

    wxButton *inputFileBrowse = new wxButton(panel, BUTTON_Browse, _T("Browse"));

    hboxInput->Add(new wxStaticText(
        panel, -1, "Input File:"));
    hboxInput->Add(inputFile);
    hboxInput->Add(inputFileBrowse);

    sizer->Add(hboxInput);


    wxBoxSizer *hboxNeighbor = new wxBoxSizer(wxHORIZONTAL);


    wxTextCtrl *neighbor = new wxTextCtrl(panel, -1, "1.0e-5");

    hboxNeighbor->Add(new wxStaticText(
        panel, -1, "Neighbor Skin Distance (m):"));
    hboxNeighbor->Add(neighbor);

    sizer->Add(hboxNeighbor);


    wxBoxSizer *hboxTimestep = new wxBoxSizer(wxHORIZONTAL);


    wxTextCtrl *timestep = new wxTextCtrl(panel, -1, "1.0");

    hboxTimestep->Add(new wxStaticText(
        panel, -1, "Time Step (s):"));
    hboxTimestep->Add(timestep);

    sizer->Add(hboxTimestep);

    sizer->Add(new wxStaticText(
        panel, -1, "", wxDefaultPosition));


    wxBoxSizer *vbox5 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox6 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox7 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox8 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox9 = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer *hboxRates = new wxBoxSizer(wxHORIZONTAL);

    vbox5->Add(new wxStaticText(
        panel, -1, ""));
    vbox5->Add(new wxStaticText(
        panel, -1, "Maximum Growth Rate (1/s):", wxDefaultPosition, wxSize(-1, 34)));
    vbox5->Add(new wxStaticText(
        panel, -1, "Decay Rate (1/s):", wxDefaultPosition, wxSize(-1, 34)));
    vbox5->Add(new wxStaticText(
        panel, -1, "Yield Coefficient (-):", wxDefaultPosition, wxSize(-1, 34)));


    wxTextCtrl *hetGrowthRate = new wxTextCtrl(panel, -1, "6.944e-5");

    wxTextCtrl *hetDecayRate = new wxTextCtrl(panel, -1, "4.6296e-6");

    wxTextCtrl *hetYieldRate = new wxTextCtrl(panel, -1, "0.61");

    vbox6->Add(new wxStaticText(
        panel, -1, "HET"), 0, wxALIGN_CENTER);
    vbox6->Add(hetGrowthRate, 0, wxALIGN_CENTER);
    vbox6->Add(hetDecayRate, 0, wxALIGN_CENTER);
    vbox6->Add(hetYieldRate, 0, wxALIGN_CENTER);


    wxTextCtrl *aobGrowthRate = new wxTextCtrl(panel, -1, "3.4722e-5");

    wxTextCtrl *aobDecayRate = new wxTextCtrl(panel, -1, "1.2731e-6");

    wxTextCtrl *aobYieldRate = new wxTextCtrl(panel, -1, "0.33");


    vbox7->Add(new wxStaticText(
        panel, -1, "AOB"), 0, wxALIGN_CENTER);
    vbox7->Add(aobGrowthRate, 0, wxALIGN_CENTER);
    vbox7->Add(aobDecayRate, 0, wxALIGN_CENTER);
    vbox7->Add(aobYieldRate, 0, wxALIGN_CENTER);


    wxTextCtrl *nobGrowthRate = new wxTextCtrl(panel, -1, "3.4722e-5");

    wxTextCtrl *nobDecayRate = new wxTextCtrl(panel, -1, "1.2731e-6");

    wxTextCtrl *nobYieldRate = new wxTextCtrl(panel, -1, "0.083");


    vbox8->Add(new wxStaticText(
        panel, -1, "NOB"), 0, wxALIGN_CENTER);
    vbox8->Add(nobGrowthRate, 0, wxALIGN_CENTER);
    vbox8->Add(nobDecayRate, 0, wxALIGN_CENTER);
    vbox8->Add(nobYieldRate, 0, wxALIGN_CENTER);


    wxTextCtrl *epsDecayRate = new wxTextCtrl(panel, -1, "1.9676e-6");

    wxTextCtrl *epsYieldRate = new wxTextCtrl(panel, -1, "0.18");


    vbox9->Add(new wxStaticText(
        panel, -1, "EPS"), 0, wxALIGN_CENTER);
    vbox9->Add(new wxStaticText(
        panel, -1, "", wxDefaultPosition, wxSize(-1, 34)), 0, wxALIGN_CENTER);
    vbox9->Add(epsDecayRate, 0, wxALIGN_CENTER);
    vbox9->Add(epsYieldRate, 0, wxALIGN_CENTER);

    hboxRates->Add(vbox5);
    hboxRates->Add(vbox6);
    hboxRates->Add(vbox7);
    hboxRates->Add(vbox8);
    hboxRates->Add(vbox9);


    sizer->Add(hboxRates);

    sizer->Add(new wxStaticText(
        panel, -1, "", wxDefaultPosition));


    wxBoxSizer *vbox10 = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer *vbox11 = new wxBoxSizer(wxVERTICAL);

    wxBoxSizer *hboxAffinities = new wxBoxSizer(wxHORIZONTAL);

    vbox10->Add(new wxStaticText(
        panel, -1, "Dissolved Oxygen Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(
        panel, -1, "Nitrite Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(
        panel, -1, "HET Carbon Source Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(
        panel, -1, "HET Nitrate Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(
        panel, -1, "AOB Ammonia Affinity (kg/m^3):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(
        panel, -1, "HET Reduction Factor (-):", wxDefaultPosition, wxSize(-1, 34)));
    vbox10->Add(new wxStaticText(
        panel, -1, "Death Factor (-):", wxDefaultPosition, wxSize(-1, 34)));


    wxBoxSizer *hboxOxygen = new wxBoxSizer(wxHORIZONTAL);

    wxTextCtrl *hetOxygen = new wxTextCtrl(panel, -1, "8.1e-4");
    wxTextCtrl *aobOxygen = new wxTextCtrl(panel, -1, "5e-4");
    wxTextCtrl *nobOxygen = new wxTextCtrl(panel, -1, "6.8e-4");

    hboxOxygen->Add(new wxStaticText(
        panel, -1, "HET:"));
    hboxOxygen->Add(hetOxygen);
    hboxOxygen->Add(new wxStaticText(
        panel, -1, "AOB:"));
    hboxOxygen->Add(aobOxygen);
    hboxOxygen->Add(new wxStaticText(
        panel, -1, "NOB:"));
    hboxOxygen->Add(nobOxygen);


    wxBoxSizer *hboxNitrite = new wxBoxSizer(wxHORIZONTAL);

    wxTextCtrl *hetNitrite = new wxTextCtrl(panel, -1, "3e-4");
    wxTextCtrl *nobNitrite = new wxTextCtrl(panel, -1, "1.3e-3");

    hboxNitrite->Add(new wxStaticText(
        panel, -1, "HET:"));
    hboxNitrite->Add(hetNitrite);
    hboxNitrite->Add(new wxStaticText(
        panel, -1, "NOB:"));
    hboxNitrite->Add(nobNitrite);

    wxTextCtrl *hetCarbon = new wxTextCtrl(panel, -1, "1e-2");
    wxTextCtrl *hetNitrate = new wxTextCtrl(panel, -1, "3e-4");
    wxTextCtrl *aobAmmonia = new wxTextCtrl(panel, -1, "1e-3");
    wxTextCtrl *hetReduction = new wxTextCtrl(panel, -1, "0.6");
    wxTextCtrl *deathFactor = new wxTextCtrl(panel, -1, "2.0");


    vbox11->Add(hboxOxygen);
    vbox11->Add(hboxNitrite);
    vbox11->Add(hetCarbon);
    vbox11->Add(hetNitrate);
    vbox11->Add(aobAmmonia);
    vbox11->Add(hetReduction);
    vbox11->Add(deathFactor);

    hboxAffinities->Add(vbox10);
    hboxAffinities->Add(vbox11);

    sizer->Add(hboxAffinities);


    wxBoxSizer *hboxEPS = new wxBoxSizer(wxHORIZONTAL);

    wxTextCtrl *epsDensity = new wxTextCtrl(panel, -1, "30");
    wxTextCtrl *epsRatio = new wxTextCtrl(panel, -1, "1.25");

    hboxEPS->Add(new wxStaticText(
        panel, -1, "EPS:"));
    hboxEPS->Add(new wxStaticText(
        panel, -1, "Density (kg/m^3):"));
    hboxEPS->Add(epsDensity);
    hboxEPS->Add(new wxStaticText(
        panel, -1, "Extraction Ratio (-):"));
    hboxEPS->Add(epsRatio);


    sizer->Add(hboxEPS);


    main->Add(sizer, 1, wxEXPAND | wxALL, 20);

    panel->SetSizer(main);

    Centre();







//     wxGrid *grid = new wxGrid( this,
//                     -1,
//                     wxPoint( 0, 240 ),
//                     wxSize( 400, 300 ) );

//     grid->CreateGrid( 100, 10 );
// // We can set the sizes of individual rows and columns
// // in pixels
// grid->SetRowSize( 0, 60 );
// grid->SetColSize( 0, 120 );
// // And set grid cell contents as strings
// grid->SetCellValue( 0, 0, "wxGrid is good" );
// // We can specify that some cells are read->only
// grid->SetCellValue( 0, 3, "This is read->only" );
// grid->SetReadOnly( 0, 3 );
// // Colours can be specified for grid cell contents
// grid->SetCellValue(3, 3, "green on grey");
// grid->SetCellTextColour(3, 3, *wxGREEN);
// grid->SetCellBackgroundColour(3, 3, *wxLIGHT_GREY);
// // We can specify the some cells will store numeric
// // values rather than strings. Here we set grid column 5
// // to hold floating point values displayed with width of 6
// // and precision of 2
// grid->SetColFormatFloat(5, 6, 2);
// grid->SetCellValue(0, 6, "3.1415");





}


// event handlers

void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
	delete inputFile;
    // true is to force the frame to close
    Close(true);
}

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
    wxMessageBox(wxString::Format
                 (
                    "Welcome to %s!\n"
                    "\n"
                    "This is the minimal wxWidgets sample\n"
                    "running under %s.",
                    wxVERSION_STRING,
                    wxGetOsDescription()
                 ),
                 "About wxWidgets minimal sample",
                 wxOK | wxICON_INFORMATION,
                 this);
}

void MyFrame::Browse(wxCommandEvent& WXUNUSED(event))
{
    wxFileDialog* OpenDialog = new wxFileDialog(
		this, "Choose a file to open", wxEmptyString, wxEmptyString, 
		"Input Files (*.in)|*.in|Text Files (*.txt)|*.txt",
		wxFD_OPEN, wxDefaultPosition);
 
	// Creates a "open file" dialog with 4 file types
	if (OpenDialog->ShowModal() == wxID_OK) // if the user click "Open" instead of "Cancel"
	{
		inputFile->SetValue(OpenDialog->GetPath());
	}
 
	// Clean up after ourselves
	OpenDialog->Destroy();
}
