import { ISample } from "../../models/view-samples.model";

export interface ILoadAllSamplesResponse {
  isSuccess: boolean;
  sampleList?: ISample[] | null;
}
