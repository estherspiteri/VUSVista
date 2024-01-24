import { IVus } from "../../models/view-vus.model";

export interface IStoreAndVerifyVusFileRequest {
  vusFile: File;
}

export interface IStoreAndVerifyVusFileResponse {
  isSuccess: boolean;
  areRsidsRetrieved: boolean;
  isClinvarAccessed: boolean;
  vusList: IVus[];
}

export interface ILoadAllVusResponse {
  isSuccess: boolean;
  vusList?: IVus[] | null;
}
