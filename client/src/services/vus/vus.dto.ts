import { IVus } from "../../models/view-vus.model";
import {
  ISamplePhenotypeSelected,
  IVusGene,
  IVusGeneSelected,
} from "../../models/vus-file-upload.model";

export interface IStoreAndVerifyVusFileRequest {
  vusFile: File;
  multipleGenesSelection?: IVusGeneSelected[];
  samplePhenotypeSelection?: ISamplePhenotypeSelected[];
}

export interface IStoreAndVerifyVusFileResponse {
  isSuccess: boolean;
  areRsidsRetrieved: boolean;
  isClinvarAccessed: boolean;
  vusList: IVus[];
  multipleGenes?: IVusGene[] | null;
}

export interface ILoadAllVusResponse {
  isSuccess: boolean;
  vusList?: IVus[] | null;
}

