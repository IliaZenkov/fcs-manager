import numpy as np
import pandas as pd


#######################################################################################################################
### FCS file standard references:    [1] https://www.bioconductor.org/packages/release/bioc/vignettes/flowCore/inst/doc/fcs3.html
###                                  [2] https://en.wikipedia.org/wiki/Flow_Cytometry_Standard
#######################################################################################################################


class convertFCS:
    """
    Contains methods and attributes to read and hold data from an FCS file

    Params:
    fcs: string: Path to FCS file
    An FCS file may be a wrapper around multiple other FCS files in an experiment with multiple samples
    sample_number: int: Specifies the current FCS file/experiment sample being converted from a multi-sample experiment

    Vars:
    self.text_keywords: Dict: keys are parsed from $TEXT fields of FCS file and hold the data under each $TEXT field

    self.data: NumPy ndarray[total number of events, number of channels/parameters]: Contains data from $BEGINDATA to $ENDDATA section of FCS file
    self.channel_names: Contains the names of the scattering and fluorescence channels included in the experiment run
    """

    def __init__(self, fcs, sample_number):
        self.data = None
        # Default FCS standard uses '$PnN' to refer to (short) names of acquisition channels
        self.channel_name = "$PnN"

        # Initialize variables to hold attributes from $TEXT fields from the FCS file
        self.data_start = -1
        self.data_end = -1
        self.channel_names = []
        self.channel_nums = []
        self.text_keywords = {}

        self.fcs = fcs
        with open(fcs, "rb") as fcs:
            fcs.seek(0, 2)
            fcs.seek(0)
            data_segments = 0
            # seek the correct data set in fcs
            data_offset = 0
            while data_segments <= sample_number:
                self.read_header(fcs, data_offset)
                self.read_text(fcs)
                data_segments += 1
                data_offset = self.text_keywords["$NEXTDATA"]
                fcs.seek(data_offset)
                self.read_data(fcs)

    def read_header(self, fcs, data_offset=0):
        """
        Records byte offsets locating data segments from the FCS file's HEADER segment
        HEADER: Dict: stores byte offsets for beginning and end of $TEXT and $DATA segments

        Relevant $TEXT segments contained in HEADER:
        $BEGINDATA Byte-offset to the beginning of the DATA segment.
        $BEGINSTEXT Byte-offset to the beginning of a supplemental TEXT segment.
        $ENDDATA Byte-offset to the end of the DATA segment.
        $ENDSTEXT Byte-offset to the end of a supplemental TEXT segment.

        :param fcs: current byte offset position in FCS file
        :param data_offset: byte offset for the HEADER segment itself
        :return: None
        """
        # Ignore first 10 bytes of HEADER contain FCS file format followed by 4 spaces
        fcs.read(10)

        for text in (
            "$BEGINSTEXT",
            "$ENDSTEXT",
            "$BEGINDATA",
            "$ENDDATA",
        ):
            text_offset = int(fcs.read(8))
            self.text_keywords[text] = text_offset + data_offset

        self.data_start = self.text_keywords["$BEGINDATA"]
        self.data_end = self.text_keywords["$BEGINDATA"]

    def read_text(self, fcs):
        """
        Extracts data from $BEGINSTEXT to $ENDSTEXT section of FCS file
        :param fcs: current byte-offset position in FCS file
        :return: None
        """
        fcs.seek(self.text_keywords["$BEGINSTEXT"], 0)
        text = fcs.read(
            self.text_keywords["$ENDSTEXT"] - self.text_keywords["$BEGINSTEXT"] + 1
        )
        text = text.decode(
            encoding="utf-8"
        )  # $TEXT keywords in the FCS standard are UTF-8 encoded

        # delimiters are sometimes used as escape chars in FCS files
        # To avoid misinterpreting escaped delimiters as actual delimiters, we should split the text on double (escaped) delimiters first, followed by splitting on actual delimiters
        # We record the data between the delimiters at either end of each text segment (text[1:-1])
        delim = text[
            0
        ]  # text segments start with the delimiter used for the whole FCS file
        text_segment_sublists = [
            text.split(delim) for text in text[1:-1].split(delim * 2)
        ]

        # Start the list of text segments with the first sublist, then extend it for each successive sublist
        text_segments = text_segment_sublists[0]
        # Now combine adjacent sublists (produced by double delimiters) back into whole text segments
        for segment_sublist in text_segment_sublists[1:]:
            text_segments[-1] += delim + segment_sublist[0]
            text_segments.extend(segment_sublist[1:])
        # In the now-flattened list text_keywords we have even elements as keys, and odd elements as values
        keys, values = text_segments[0::2], text_segments[1::2]
        # Finally, convert raw text data to a dictionary and add it to our list of $TEXT data in self.text_keywords
        self.text_keywords.update(dict(zip(keys, values)))

        # Don't forget to convert keys in self.text_keywords which encode bit values into integers
        # $PnB is FCS standard for bits reserved for each parameter (channel)
        bit_keys = [f"$P{i}B" for i in self.channel_nums]
        bit_keys.extend(
            ["$NEXTDATA", "$PAR", "$TOT"]
        )  # These text keywords giving byte offsets are also encoded as bits
        for key in bit_keys:
            value = self.text_keywords[key]
            self.text_keywords[key] = int(value)

        # Find how many flow acquisition channels were enabled (i.e. included parameters as in FCS standard)
        # $PAR: FCS standard keyword for parameters (channels) recorded for each event i.e. included channels numbers
        # Channel numbers start from 1 on the FACSAria, FACSCalibur, and Guava easyCyte flow cytometers (may not be true for other cytometer models/setups!)
        num_channels = int(self.text_keywords["$PAR"])
        self.channel_nums = range(1, num_channels + 1)

        # Find names of acquisition channels from $PnN keyword (i.e. parameter names as in FCS standard)
        self.channel_names = [self.text_keywords[f"$P{i}N"] for i in self.channel_nums]

        # Don't forget to convert keys in self.text_keywords which encode bit values into integers
        # $PnB is FCS standard for bits reserved for each parameter (channel)
        bit_keys = [f"$P{i}B" for i in self.channel_nums]
        bit_keys.extend(
            ["$NEXTDATA", "$PAR", "$TOT"]
        )  # These text keywords giving byte offsets are also encoded as bits
        for key in bit_keys:
            value = self.text_keywords[key]
            self.text_keywords[key] = int(value)

        ##### Update $BEGINDATA segments if needed
        if self.data_start == 0:
            self.data_start = int(text["$BEGINDATA"])
        if self.data_end == 0:
            self.data_end = int(text["$ENDDATA"])

    def read_data(self, fcs):
        """
        Extracts data from $BEGINDATA to $ENDDATA section of FCS file and records to self.data 
        :param fcs:
        :return:
        """
        text = self.text_keywords

        # $TOT: FCS keyword for number of events in the data set
        num_events = text["$TOT"]
        # $PAR: FCS keyword for Parameters (channels) recorded for each event i.e. enabled channels in the experiment
        num_params = text["$PAR"]

        # The following block attempts to catch an FCS file with mixed datatypes for different acquisition channels
        # This should not be the case, but if it is - then this will prevent a corrupt output from a valid FCS file making its way into analysis
        # $PnB: FCS standard for number of ****BITS**** reserved for parameter number n i.e. for data recorded by each channel
        reserved_bytes = [int(text[f"$P{i}B"] / 8) for i in self.channel_nums]
        # FCS $DATATYPE keyword uses $DATATYPE='D' or $DATATYPE='F' to refer to floats, and $DATATYPE='I' to integers
        # We can use a dictionary in a clever way to convert the FCS datatypes to NumPy compatible dtypes, like so:
        param_dtype = {"I": "u", "D": "f", "F": "f"}[text["$DATATYPE"]]
        parameter_data_dtypes = [
            f"{param_dtype}{parameter_bytes}" for parameter_bytes in reserved_bytes
        ]
        if len(set(parameter_data_dtypes)) > 1:
            raise NameError(
                "Mixed data types across FCS channels; check instrument settings"
            )

        fcs.seek(self.data_start, 0)  # Go back to $BEGINDATA

        # Datatypes of all channels should be the same, so pick dtype from first channel
        dtype = parameter_data_dtypes[0]
        data = np.fromfile(fcs, dtype=dtype, count=num_events * num_params)
        data = data.reshape(num_events, num_params)
        self.data = data


def convertDF(fcs, sample_number=0):
    """
    Runs convertFCS to create an FCSdata object containing flow data from FCS file
    Converts data values and keys (channel names) to a PANDAS dataframe
    :param sample_number: number of samples represented in the FCS file 
    :return: PANDAS dataframe object
    """
    # Create an FCSdata object holding data from the FCS file
    FCSdata = convertFCS(fcs, sample_number)
    # Convert the cleaned FCS data into a PANDAS dataframe
    # Clean up the df by converting to int since values after decimal won't change the analysis
    df = pd.DataFrame(FCSdata.data, columns=FCSdata.channel_names).astype(int)
    return df
